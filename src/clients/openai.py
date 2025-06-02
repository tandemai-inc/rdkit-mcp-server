import asyncio
import json
import logging

from agents import Agent, Runner, gen_trace_id, trace
from agents.mcp import MCPServer, MCPServerSse
from agents.model_settings import ModelSettings

logger = logging.getLogger(__name__)

MCP_URL = "http://localhost:8000/sse"
MCP_NAME = "RDKIT MCP Server"
OPENAI_TRACE_URL = "https://platform.openai.com/traces/trace?trace_id={}"

AGENT_INSTRUCTIONS = (
    "You are an agent that aids scientists working in the field of chemistry. "
    "You have access to a set of tools that can be used to perform various cheminformatics tasks. "
    "Use the tools to answer the user's questions. All numeric values in response must be based on the output of the tools. "
    "The final output will be read in a terminal; do not use Markdown or any other formatting. "
)

# Default prompt makes testing more convenient
DEFAULT_PROMPT = 'What is the molecular weight of `CC(=O)NC1=CC=C(C=C1)O`'


async def run(mcp_server: MCPServer, prompt: str = None) -> Runner:
    prompt = prompt or ""
    agent = Agent(
        name="RDKIT Agent",
        instructions=AGENT_INSTRUCTIONS,
        mcp_servers=[mcp_server],
        model_settings=ModelSettings(tool_choice="required"),
    )
    result: Runner = await Runner.run(starting_agent=agent, input=prompt)
    # output = format_final_output(result)
    return result


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
                result: Runner = await run(server, prompt)
                print(result.final_output)


def parse_function_calls(runner: Runner):
    """Parse function calls and outputs from the runner and return as a list."""
    input_list = runner.to_input_list()
    function_calls = []
    for input_item in input_list:
        function_types = ['function_call', 'function_call_output']
        if 'type' in input_item and input_item['type'] in function_types:
            function_calls.append(input_item)
    return function_calls


def format_final_output(runner: Runner) -> str:
    """Format the final output of the runner in a human readable format."""
    final_output = f'FINAL OUTPUT: {runner.final_output}\n\n'
    function_calls = parse_function_calls(runner)
    for call in function_calls:
        if 'arguments' in call:
            final_output += f"Function Call: {call.get('name', 'unknown')}\n"
            try:
                args = json.loads(call['arguments'])
                arg_string = '\n'.join(f"{k}: {v}" for k, v in args.items())
            except (json.JSONDecodeError, TypeError):
                arg_string = str(call['arguments'])
            final_output += f"Arguments: {arg_string}\n"
        if 'output' in call:
            try:
                output = json.loads(call['output'])
                output_str = '\n'.join(f"{k}: {v}" for k, v in output.items())
            except (json.JSONDecodeError, TypeError):
                output_str = str(call['output'])
            final_output += f"Output: {output_str}\n"
        final_output += "\n"
    return final_output

if __name__ == "__main__":
    asyncio.run(main())
