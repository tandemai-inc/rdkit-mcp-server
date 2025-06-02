from typing import Dict, Any, Optional, Union, List
from src.clients.openai import MCP_URL, MCP_NAME, OPENAI_TRACE_URL, run

from agents import gen_trace_id, trace
from agents.mcp import MCPServerSse


async def call_llm(prompt: str) -> str:
    async with MCPServerSse(
            name=MCP_NAME,
            params={
                "url": MCP_URL,
            },
        ) as server:
            trace_id = gen_trace_id()
            with trace(workflow_name=prompt, trace_id=trace_id):
                print(f"View trace: {OPENAI_TRACE_URL.format(trace_id)}\n")
                final_output = await run(server, prompt)
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
    # Note: The prompt may be in JSON format, so you might need to parse it.
    # For example, if the prompt is a JSON string representing a conversation:
    # prompt = '[{"role": "user", "content": "Hello, world!"}]'
    # You would parse it like this:
    # prompt = json.loads(prompt)

    # The 'options' parameter contains additional configuration for the API call.
    config = options.get('config', None)
    # additional_option = config.get('additionalOption', None)

    # The 'context' parameter provides info about which vars were used to create the final prompt.
    # user_variable = context['vars'].get('userVariable', None)

    # The prompt is the final prompt string after the variables have been processed.
    # Custom logic to process the prompt goes here.
    # For instance, you might call an external API or run some computations.
    # TODO: Replace with actual LLM API implementation.
    output = await call_llm(prompt)

    # The result should be a dictionary with at least an 'output' field.
    result = {
        "output": output,
    }
    # if some_error_condition:
    #     result['error'] = "An error occurred during processing"

    # if token_usage_calculated:
    #     # If you want to report token usage, you can set the 'tokenUsage' field.
    #     result['tokenUsage'] = {
    #         "total": token_count, "prompt": prompt_token_count, "completion": completion_token_count}

    # if failed_guardrails:
    #     # If guardrails triggered, you can set the 'guardrails' field.
    #     result['guardrails'] = {"flagged": True}

    return result


def call_embedding_api(prompt: str) -> ProviderEmbeddingResponse:
    # Returns ProviderEmbeddingResponse
    pass


def call_classification_api(prompt: str) -> ProviderClassificationResponse:
    # Returns ProviderClassificationResponse
    pass
