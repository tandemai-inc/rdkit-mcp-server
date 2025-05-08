import asyncio

from agents import Agent, Runner, gen_trace_id, trace
from agents.mcp import MCPServer, MCPServerSse
from agents.model_settings import ModelSettings


async def run(mcp_server: MCPServer, prompt: str = None):
    prompt = prompt or ""
    agent = Agent(
        name="Assistant",
        instructions="Use the tools to answer the questions.",
        mcp_servers=[mcp_server],
        model_settings=ModelSettings(tool_choice="required"),
    )
    print(f"Running: {prompt}")
    result = await Runner.run(starting_agent=agent, input=prompt)
    print(result.final_output)


async def main():
    # Create a while loop that requests a prompt from a user, and then sends it to the agent to process
    while True:
        prompt = input("Enter a prompt (or 'exit' to quit): ")
        if prompt.lower() == "exit":
            break
        # Create a new MCPServer instance
        async with MCPServerSse(
            name="SSE Python Server",
            params={
                "url": "http://localhost:8000/sse",
            },
        ) as server:
            trace_id = gen_trace_id()
            with trace(workflow_name="MCP Example", trace_id=trace_id):
                print(f"View trace: https://platform.openai.com/traces/trace?trace_id={trace_id}\n")
                await run(server, prompt)

if __name__ == "__main__":
    asyncio.run(main())
