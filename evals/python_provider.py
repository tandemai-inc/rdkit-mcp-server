import json
import logging
import openai
import time

from agents import Runner, gen_trace_id, trace
from agents.mcp import MCPServerSse
from typing import Dict, Any, Optional, Union, List
from rdkit_mcp_clients.openai import MCP_URL, MCP_NAME, create_agent


logger = logging.getLogger(__name__)


async def call_llm(prompt: str, model=None, use_mcp=True) -> Runner:
    async with MCPServerSse(
        name=MCP_NAME,
        params={
            "url": MCP_URL,
        },
    ) as server:
        trace_id = gen_trace_id()
        with trace(workflow_name=prompt, trace_id=trace_id):
            kwargs = {}
            if use_mcp:
                kwargs["mcp_server"] = server
            agent = create_agent(model=model, **kwargs)
            result: Runner = None
            retry_attempts = 5
            for i in range(retry_attempts):
                try:
                    result: Runner = await Runner.run(starting_agent=agent, input=prompt)
                    break
                except openai.RateLimitError:
                    wait = 2 ** (i + 5)  # Exponential backoff. Start at 32 seconds
                    logger.error(f"Rate limit hit. Retrying in {wait} seconds...")
                    time.sleep(wait)
            return result


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
                if isinstance(args, list):
                    arg_string = '\n'.join({f"arg_{i}": v for i, v in enumerate(args)})
                else:
                    arg_string = '\n'.join(f"{k}: {v}" for k, v in args.items())
            except (json.JSONDecodeError, TypeError):
                arg_string = str(call['arguments'])
            final_output += f"Arguments: {arg_string}\n"
        if 'output' in call:
            try:
                output = json.loads(call['output'])
                if isinstance(output, dict):
                    output_str = '\n'.join(f"{k}: {v}" for k, v in output.items())
                else:
                    output_str = str(output)
            except (json.JSONDecodeError, TypeError):
                output_str = str(call['output'])
            final_output += f"Output: {output_str}\n"
        final_output += "\n"
    return final_output


class ProviderOptions:
    id: Optional[str]
    config: Optional[Dict[str, Any]]


class CallApiContextParams:
    vars: Dict[str, str]


class TokenUsage:
    total: int
    prompt: int
    completion: int


class ProviderResponse:
    output: Optional[Union[str, Dict[str, Any]]]
    error: Optional[str]
    tokenUsage: Optional[TokenUsage]
    cost: Optional[float]
    cached: Optional[bool]
    logProbs: Optional[List[float]]


class ProviderEmbeddingResponse:
    embedding: List[float]
    tokenUsage: Optional[TokenUsage]
    cached: Optional[bool]


class ProviderClassificationResponse:
    classification: Dict[str, Any]
    tokenUsage: Optional[TokenUsage]
    cached: Optional[bool]


async def call_api(prompt: str, options: Dict[str, Any], context: Dict[str, Any]) -> ProviderResponse:
    # The 'options' parameter contains additional configuration for the API call.
    config = options.get('config', None)
    model = config.get('model', None)
    use_mcp = config.get('use_mcp', True)

    result: Runner = await call_llm(prompt, model=model, use_mcp=use_mcp)
    if result:
        final_output = format_final_output(result)
    else:
        final_output = "Error: No result returned from the LLM."
    response = {
        "output": final_output
    }
    return response


def call_embedding_api(prompt: str) -> ProviderEmbeddingResponse:
    # Returns ProviderEmbeddingResponse
    pass


def call_classification_api(prompt: str) -> ProviderClassificationResponse:
    # Returns ProviderClassificationResponse
    pass
