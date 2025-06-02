import logging
from agents import Runner, gen_trace_id, trace
from agents.mcp import MCPServerSse
from typing import Dict, Any, Optional, Union, List
from src.clients.openai import MCP_URL, MCP_NAME, OPENAI_TRACE_URL, run, format_final_output
logger = logging.getLogger(__name__)


async def call_llm(prompt: str, model=None) -> Runner:
    async with MCPServerSse(
        name=MCP_NAME,
        params={
            "url": MCP_URL,
        },
    ) as server:
        trace_id = gen_trace_id()
        with trace(workflow_name=prompt, trace_id=trace_id):
            print(f"View trace: {OPENAI_TRACE_URL.format(trace_id)}\n")
            result: Runner = await run(server, prompt, model=model)
    return result


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

    result: Runner = await call_llm(prompt, model=model)
    final_output = format_final_output(result)
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
