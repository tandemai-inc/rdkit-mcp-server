"""Task function for RDKit MCP evaluations using pydantic-ai."""

import asyncio
import logging
from dataclasses import dataclass, field
from typing import Any

from pydantic_ai import Agent
from pydantic_ai.mcp import MCPServerSSE
from pydantic_ai.messages import ModelResponse
from pydantic_ai.exceptions import UnexpectedModelBehavior

logger = logging.getLogger(__name__)

MCP_URL = "http://localhost:8000/sse"

AGENT_INSTRUCTIONS = (
    "You are an agent that aids scientists working in the field of chemistry. "
    "You have access to a set of tools that can be used to perform various cheminformatics tasks. "
    "Use the tools to answer the user's questions. All numeric values in response must be based on the output of the tools. "
    "The final output will be read in a terminal; do not use Markdown or any other formatting. "
    "If the final output is a file, use the write_file tool to write the file and return the file path in the final output."
)


@dataclass
class ToolCall:
    """A single tool call made during task execution."""

    tool_name: str
    args: dict[str, Any]
    tool_call_id: str


@dataclass
class TaskInput:
    """Input for the evaluation task."""

    prompt: str
    model: str = "openai:gpt-4o"
    use_mcp: bool = True


@dataclass
class TaskOutput:
    """Output from the evaluation task, including tool call history."""

    text: str
    tool_calls: list[ToolCall] = field(default_factory=list)

    def __str__(self) -> str:
        """Return the text output for string representation."""
        return self.text


def _extract_tool_calls(messages: list) -> list[ToolCall]:
    """Extract all tool calls from the message history."""
    tool_calls = []
    for msg in messages:
        if isinstance(msg, ModelResponse):
            for part in msg.parts:
                if part.part_kind == "tool-call":
                    tool_calls.append(
                        ToolCall(
                            tool_name=part.tool_name,
                            args=part.args if isinstance(part.args, dict) else {},
                            tool_call_id=part.tool_call_id,
                        )
                    )
    return tool_calls


async def run_task_async(inputs: TaskInput) -> TaskOutput:
    """Execute a prompt against the RDKit MCP server via pydantic-ai Agent."""
    server = MCPServerSSE(MCP_URL)

    if inputs.use_mcp:
        agent = Agent(
            inputs.model,
            system_prompt=AGENT_INSTRUCTIONS,
            toolsets=[server],
        )
    else:
        agent = Agent(
            inputs.model,
            system_prompt=AGENT_INSTRUCTIONS,
        )

    retry_attempts = 5
    last_error: Exception | None = None

    async with agent:
        for i in range(retry_attempts):
            try:
                result = await agent.run(inputs.prompt)
                tool_calls = _extract_tool_calls(result.all_messages())
                return TaskOutput(text=str(result.output), tool_calls=tool_calls)
            except UnexpectedModelBehavior as e:
                last_error = e
                wait = 2 ** (i + 5)  # Exponential backoff starting at 32s
                logger.error(f"Model error. Retrying in {wait} seconds...")
                await asyncio.sleep(wait)
            except Exception as e:
                # Check if it's a rate limit error
                if "rate" in str(e).lower() or "429" in str(e):
                    last_error = e
                    wait = 2 ** (i + 5)
                    logger.error(f"Rate limit hit. Retrying in {wait} seconds...")
                    await asyncio.sleep(wait)
                else:
                    raise

    if last_error:
        return TaskOutput(text=f"Error after {retry_attempts} retries: {last_error}")
    return TaskOutput(text="Error: No result returned from the agent.")


def run_task_sync(inputs: TaskInput) -> TaskOutput:
    """Synchronous wrapper for evaluate_sync compatibility."""
    return asyncio.run(run_task_async(inputs))
