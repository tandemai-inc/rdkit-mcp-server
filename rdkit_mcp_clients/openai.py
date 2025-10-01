import asyncio
import logging
import os
from typing import Annotated
import openai
import time

from agents import Agent, Runner, gen_trace_id, trace
from agents.mcp import MCPServer, MCPServerSse
from agents.model_settings import ModelSettings

from agents.tool import function_tool
import base64

from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)

MCP_URL = "http://localhost:8000/sse"
MCP_NAME = "RDKIT MCP Server"
OPENAI_TRACE_URL = "https://platform.openai.com/traces/trace?trace_id={}"

AGENT_INSTRUCTIONS = (
    "You are an agent that aids scientists working in the field of chemistry. "
    "You have access to a set of tools that can be used to perform various cheminformatics tasks. "
    "Use the tools to answer the user's questions. All numeric values in response must be based on the output of the tools. "
    "The final output will be read in a terminal; do not use Markdown or any other formatting. "
    "If the final output is a file, use the write_file tool to write the file and return the file path in the final output. "
)


class EncodedFileModel(BaseModel):
    filename: Annotated[str, Field(description="Name of the file.")]
    content: Annotated[str, Field(description="Base 64 encoded bytes containing contents of the file.")]


@function_tool(description_override="Reads the contents of a file at the given absolute file path.")
def read_file(filepath: str) -> str:
    with open(filepath, "r") as f:
        return f.read()


@function_tool(
    description_override=(
        "Writes a base64-encoded file to the specified directory. "
        "Provide the target directory and an EncodedFileModel containing the filename and base64-encoded content. "
        "Returns the absolute path to the written file."
    )
)
def write_file(directory: str, encoded_file: EncodedFileModel) -> str:
    try:
        decoded_content = base64.b64decode(encoded_file.content)
    except Exception as e:
        raise ValueError(f"Failed to decode file content: {e}")

    filepath = os.path.join(directory, encoded_file.filename)
    with open(filepath, "wb") as f:
        f.write(decoded_content)
    return f"Successfully wrote file to {filepath}"


def create_agent(mcp_server: MCPServer = None, model: str = None) -> Agent:
    """Create an agent with the specified MCP server and model."""
    mcp_servers = []
    if mcp_server:
        mcp_servers.append(mcp_server)

    agent = Agent(
        name="RDKIT Agent",
        instructions=AGENT_INSTRUCTIONS,
        mcp_servers=mcp_servers,
        model_settings=ModelSettings(tool_choice="auto"),
        model=model
    )
    agent.tools.append(read_file)
    agent.tools.append(write_file)
    return agent


async def main(prompt: str = None, model: str = "4o-mini"):
    prompt = prompt or ""

    async with MCPServerSse(
        name=MCP_NAME,
        params={"url": MCP_URL},
    ) as server:
        agent = create_agent(mcp_server=server, model=model)

        conversation_history = []
        while True:
            if not prompt:
                prompt = input("Enter a prompt or 'exit': ")
            if prompt.lower().strip() == "exit":
                break
            # Append user's message to history
            conversation_history.append({"role": "user", "content": prompt})

            # Create a single input string from history
            full_prompt = "\n".join(
                [f"{msg['role']}: {msg['content']}" for msg in conversation_history]
            )

            trace_id = gen_trace_id()
            with trace(workflow_name=prompt, trace_id=trace_id):
                print(f"View trace: {OPENAI_TRACE_URL.format(trace_id)}\n")
                retry_attempts = 5
                result = None
                for i in range(retry_attempts):
                    try:
                        result: Runner = await Runner.run(starting_agent=agent, input=full_prompt)
                        break
                    except openai.RateLimitError:
                        wait = 2 ** (i + 5)
                        logger.error(f"Rate limit hit. Retrying in {wait} seconds...")
                        time.sleep(wait)

                if result:
                    print(f"\n{result.final_output}\n")
                    # Add assistant's response to history
                    conversation_history.append({"role": "assistant", "content": result.final_output})
                prompt = ""

if __name__ == "__main__":
    asyncio.run(main())
