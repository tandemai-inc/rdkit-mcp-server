"""Task function for RDKit MCP evaluations using pydantic-ai."""

import asyncio
import logging
import time
from dataclasses import dataclass
from typing import Any

from pydantic_ai import Agent
from pydantic_ai.mcp import MCPServerSSE
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
class TaskInput:
    """Input for the evaluation task."""

    prompt: str
    model: str = "openai:gpt-4o-mini"
    use_mcp: bool = True


async def run_task_async(inputs: TaskInput) -> str:
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
                return str(result.output)
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
        return f"Error after {retry_attempts} retries: {last_error}"
    return "Error: No result returned from the agent."


def run_task_sync(inputs: TaskInput) -> str:
    """Synchronous wrapper for evaluate_sync compatibility."""
    return asyncio.run(run_task_async(inputs))
