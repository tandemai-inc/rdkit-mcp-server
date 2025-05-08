import asyncio

from agents import Agent, Runner, gen_trace_id, trace
from agents.mcp import MCPServer, MCPServerSse
from agents.model_settings import ModelSettings


AGENT_INSTRUCTIONS = (
    "You are an agent that leverages the RDkit library to aid scientists in the field of chemistry. "
    "Use the tools to answer the user's questions. All numeric values in response must be based on the output of the tools. "
    "The final output will be read in a terminal; do not use Markdown or any other formatting. "
)


async def run(mcp_server: MCPServer, prompt: str = None):
    prompt = prompt or ""
    agent = Agent(
        name="RDKIT Agent",
        instructions=AGENT_INSTRUCTIONS,
        mcp_servers=[mcp_server],
        model_settings=ModelSettings(tool_choice="required"),
    )
    result = await Runner.run(starting_agent=agent, input=prompt)
    print(result.final_output)

MCP_URL = "http://localhost:8000/sse"
MCP_NAME = "RDKIT MCP Server"
OPENAI_TRACE_URL = "https://platform.openai.com/traces/trace?trace_id={}"

# Default prompt makes testing more convenient
DEFAULT_PROMPT = 'Get basic properties of SMILES CC(=O)NC1=CC=C(C=C1)O'

async def main():
    # Create a while loop that requests a prompt from a user, and then sends it to the agent to process
    while True:
        prompt = input("Enter a prompt (or 'exit' to quit): ")
        if prompt.lower() == "exit":
            break
        if not prompt:
            prompt = DEFAULT_PROMPT

        async with MCPServerSse(
            name=MCP_NAME,
            params={
                "url": MCP_URL,
            },
        ) as server:
            trace_id = gen_trace_id()
            with trace(workflow_name=prompt, trace_id=trace_id):
                print(f"View trace: {OPENAI_TRACE_URL.format(trace_id)}\n")
                await run(server, prompt)

if __name__ == "__main__":
    asyncio.run(main())
