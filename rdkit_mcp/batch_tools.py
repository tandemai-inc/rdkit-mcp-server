from typing import Any, Dict, List, Optional
from pydantic import BaseModel
import asyncio

from mcp.server.fastmcp import Context

from .decorators import rdkit_tool


class BatchItemResult(BaseModel):
    ok: bool
    input: Optional[Any] = None
    output: Optional[Any] = None
    error: Optional[str] = None


class BatchMapResponse(BaseModel):
    tool_name: str
    results: List[BatchItemResult]


@rdkit_tool()
async def batch_map(
    tool_name: str,
    inputs: List[Dict[str, Any]],
    ctx: Context,
    concurrency: int = 10,
    fail_fast: bool = False,
    include_input: bool = True,
) -> BatchMapResponse:
    """
    Run a single MCP tool over a list of inputs in batch. Use this tool when you need to
    apply the same operation to multiple molecules efficiently in a single call.

    Example: To calculate molecular weight for 3 molecules, call batch_map with:
        tool_name="MolWt"
        inputs=[{"smiles": "CCO"}, {"smiles": "CC(=O)O"}, {"smiles": "C1=CC=CC=C1"}]

    Args:
        tool_name: Name of the MCP tool to call for each input (e.g., "MolWt", "NumRotatableBonds")
        inputs: List of argument dictionaries to map over (each dict is passed as tool arguments)
        ctx: MCP context (automatically injected)
        concurrency: Max concurrent calls (1-100, default 10)
        fail_fast: Stop on first error (default False)
        include_input: Include original input in each result item (default True)

    Returns:
        BatchMapResponse with results for each input
    """
    if not 1 <= concurrency <= 100:
        raise ValueError("concurrency must be between 1 and 100")

    sem = asyncio.Semaphore(concurrency)

    async def run_one(x: Dict[str, Any]) -> BatchItemResult:
        async with sem:
            try:
                out = await ctx.fastmcp.call_tool(tool_name, x)
                return BatchItemResult(ok=True, input=x if include_input else None, output=out)
            except Exception as e:
                return BatchItemResult(ok=False, input=x if include_input else None, error=str(e))

    results: List[BatchItemResult] = []
    for coro in asyncio.as_completed([run_one(x) for x in inputs]):
        item = await coro
        results.append(item)
        if fail_fast and not item.ok:
            break

    return BatchMapResponse(tool_name=tool_name, results=results)
