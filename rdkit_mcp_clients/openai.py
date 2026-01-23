import asyncio
import logging

from pydantic_ai import Agent, RunContext
from pydantic_ai.mcp import MCPServerSSE
from pydantic_ai.exceptions import UnexpectedModelBehavior


logger = logging.getLogger(__name__)

MCP_URL = "http://localhost:8000/sse"

AGENT_INSTRUCTIONS = (
    "You are an agent that aids scientists working in the field of chemistry. "
    "You have access to a set of tools that can be used to perform various cheminformatics tasks. "
    "Use the tools to answer the user's questions. All numeric values in response must be based on the output of the tools. "
    "The final output will be read in a terminal; do not use Markdown or any other formatting. "
    "If the final output is a file, use the write_file tool to write the file and return the file path in the final output. "
)


def create_agent(mcp_url: str = None, model: str = None) -> Agent:
    """Create an agent with the specified MCP server and model."""
    toolsets = []
    if mcp_url:
        toolsets.append(MCPServerSSE(mcp_url))

    # Format model string for pydantic-ai (requires provider prefix)
    model_str = model or "gpt-4o"
    if not model_str.startswith("openai:"):
        model_str = f"openai:{model_str}"

    agent: Agent[None, str] = Agent(
        model_str,
        system_prompt=AGENT_INSTRUCTIONS,
        toolsets=toolsets,
    )

    # Register local tools
    @agent.tool
    async def read_file(ctx: RunContext[None], filepath: str) -> str:
        """Reads the contents of a file at the given absolute file path."""
        with open(filepath, "r") as f:
            return f.read()

    @agent.tool
    async def write_file(ctx: RunContext[None], filepath: str, file_content: str) -> str:
        """Writes file_contents to the specified filepath. Returns the absolute path to the written file."""
        with open(filepath, "w") as f:
            f.write(file_content)
        return filepath

    return agent


async def main(prompt: str = None, model: str = "gpt-4o-mini"):
    prompt = prompt or ""

    agent = create_agent(mcp_url=MCP_URL, model=model)

    conversation_history = []

    async with agent:
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

            retry_attempts = 5
            result = None
            last_error = None

            for i in range(retry_attempts):
                try:
                    result = await agent.run(full_prompt)
                    break
                except UnexpectedModelBehavior as e:
                    last_error = e
                    wait = 2 ** (i + 5)
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

            if result:
                output = str(result.output)
                print(f"\n{output}\n")
                # Add assistant's response to history
                conversation_history.append({"role": "assistant", "content": output})
            elif last_error:
                print(f"\nError after {retry_attempts} retries: {last_error}\n")

            prompt = ""


if __name__ == "__main__":
    asyncio.run(main())
