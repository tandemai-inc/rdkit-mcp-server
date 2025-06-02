import logging
from agents import Runner, gen_trace_id, trace
from agents.mcp import MCPServerSse
from openai.types.responses import ResponseFunctionToolCall, ResponseOutputMessage
from typing import Dict, Any, Optional, Union, List
from src.clients.openai import MCP_URL, MCP_NAME, OPENAI_TRACE_URL, run

logger = logging.getLogger(__name__)

async def call_llm(prompt: str) -> Runner:
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
    result: Runner = await call_llm(prompt)

    # Parse tool calls and return in result
    tool_calls = []
    for model_response in result.raw_responses:
        resp_outputs: List[ResponseFunctionToolCall] = model_response.output
        for tool_call in resp_outputs:
            if isinstance(tool_call, ResponseFunctionToolCall):
                tool_calls.append(tool_call.__dict__)
            elif isinstance(tool_call, ResponseOutputMessage):
                logger.debug("Received a ResponseOutputMessage, which is not a tool call.")
            else:
                logger.warning(
                    f"Unexpected output type: {type(tool_call)}. Expected ResponseFunctionToolCall."
                )

    # The result should be a dictionary with at least an 'output' field.
    response = {
        "output":{
            "final_output": result.final_output,
            "tool_calls": tool_calls,
        }
    }

    return response


def call_embedding_api(prompt: str) -> ProviderEmbeddingResponse:
    # Returns ProviderEmbeddingResponse
    pass


def call_classification_api(prompt: str) -> ProviderClassificationResponse:
    # Returns ProviderClassificationResponse
    pass
