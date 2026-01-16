from typing import Any, List, Optional
from pydantic import BaseModel, Field
import asyncio

# Assuming FastMCP-ish patterns; adapt names to your MCP lib
from mcp.server.fastmcp import FastMCP

mcp = FastMCP("chem-tools")

class BatchMapRequest(BaseModel):
    tool_name: str = Field(..., description="Name of the MCP tool to call for each input")
    inputs: List[Any] = Field(..., description="List of inputs to map over")
    concurrency: int = Field(10, ge=1, le=100, description="Max concurrent calls")
    fail_fast: bool = Field(False, description="Stop on first error")
    include_input: bool = Field(True, description="Include original input in each result item")

class BatchItemResult(BaseModel):
    ok: bool
    input: Optional[Any] = None
    output: Optional[Any] = None
    error: Optional[str] = None

class BatchMapResponse(BaseModel):
    tool_name: str
    results: List[BatchItemResult]

@mcp.tool()
async def batch_map(req: BatchMapRequest) -> BatchMapResponse:
    """
    Run a single-input MCP tool over a list of inputs.
    """
    # Resolve the tool function by name from the MCP registry.
    # How you do this depends on your MCP framework.
    tool_fn = mcp.tools.get(req.tool_name)  # <-- adapt if your registry differs
    if tool_fn is None:
        raise ValueError(f"Unknown tool: {req.tool_name}")

    sem = asyncio.Semaphore(req.concurrency)

    async def run_one(x: Any) -> BatchItemResult:
        async with sem:
            try:
                out = await tool_fn(x)  # assumes tool signature tool(input) -> output
                return BatchItemResult(ok=True, input=x if req.include_input else None, output=out)
            except Exception as e:
                return BatchItemResult(ok=False, input=x if req.include_input else None, error=str(e))

    results: List[BatchItemResult] = []
    for coro in asyncio.as_completed([run_one(x) for x in req.inputs]):
        item = await coro
        results.append(item)
        if req.fail_fast and not item.ok:
            break

    return BatchMapResponse(tool_name=req.tool_name, results=results)
