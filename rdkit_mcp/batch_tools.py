from typing import Annotated, Any, Dict, List, Optional
from pydantic import BaseModel, Field
import asyncio

from mcp.server.fastmcp import Context

from .decorators import rdkit_tool


def optimal_concurrency(n_items: int, max_concurrency: int = 50) -> int:
    """
    Calculate optimal concurrency based on dataset size.

    Rules of thumb:
    - Small datasets (<=5): run all concurrently
    - Medium datasets (6-20): cap at 10
    - Larger datasets (21-50): cap at 20
    - Large datasets (50+): cap at max_concurrency (default 50)

    Beyond ~50 concurrent tasks, scheduling overhead typically
    exceeds parallelism benefits.
    """
    if n_items <= 5:
        return n_items
    elif n_items <= 20:
        return min(10, n_items)
    elif n_items <= 50:
        return min(20, n_items)
    else:
        return min(max_concurrency, n_items)


async def resolve_tool_name(tool_name: str, ctx: Context, sample_input: Optional[Dict[str, Any]] = None) -> str:
    """
    Resolve a potentially prefixed tool name to the actual tool name on the MCP server.

    The tool_name might come with a prefix (e.g., "rdkit_smiles_to_mol") but
    the FastMCP server only knows the unprefixed name (e.g., "smiles_to_mol").
    This function handles disambiguation when multiple tools match.

    Args:
        tool_name: The tool name (possibly with prefix)
        ctx: MCP context for accessing available tools
        sample_input: Optional sample input dict to help disambiguate based on parameters

    Returns:
        The actual tool name on the MCP server

    Raises:
        ValueError: If tool is not found or ambiguous matches cannot be resolved
    """
    available_tools = await ctx.fastmcp.list_tools()
    tool_names_set = {t.name for t in available_tools}
    tools_by_name = {t.name: t for t in available_tools}

    # If exact match exists, use it
    if tool_name in tool_names_set:
        return tool_name

    # Find all tools that the provided tool_name ends with
    matching_tools = [
        known_tool for known_tool in tool_names_set
        if tool_name.endswith(known_tool)
    ]

    if len(matching_tools) == 0:
        raise ValueError(f"Tool '{tool_name}' not found. Available tools: {sorted(tool_names_set)}")
    elif len(matching_tools) == 1:
        return matching_tools[0]

    # Multiple matches - disambiguate by checking parameters
    input_params = set(sample_input.keys()) if sample_input else set()

    best_match = None
    best_match_score = -1

    for candidate_tool_name in matching_tools:
        tool_info = tools_by_name[candidate_tool_name]
        # Get required parameters from the tool's input schema
        tool_params = set()
        if hasattr(tool_info, 'inputSchema') and tool_info.inputSchema:
            schema = tool_info.inputSchema
            if isinstance(schema, dict) and 'properties' in schema:
                tool_params = set(schema['properties'].keys())

        # Score based on parameter overlap and name length
        if tool_params and input_params:
            overlap = len(input_params & tool_params)
            # Prefer longer tool names (more specific) as tiebreaker
            score = (overlap, len(candidate_tool_name))
        else:
            # No schema info or no input params, prefer longer name (more specific match)
            score = (0, len(candidate_tool_name))

        if score > best_match_score:
            best_match_score = score
            best_match = candidate_tool_name

    if best_match:
        return best_match

    raise ValueError(
        f"Ambiguous tool name '{tool_name}'. Multiple matches found: {matching_tools}. "
        f"Cannot disambiguate based on parameters."
    )


class BatchItemResult(BaseModel):
    ok: bool
    input: Optional[Any] = None
    output: Optional[Any] = None
    error: Optional[str] = None


class BatchMapResponse(BaseModel):
    tool_name: str
    concurrency: int
    results: List[BatchItemResult]


@rdkit_tool()
async def batch_map(
    tool_name: str,
    inputs: Annotated[List[Dict[str, Any]], Field(max_length=128)],
    ctx: Context,
    concurrency: Optional[int] = None,
    fail_fast: bool = False,
    include_input: bool = True,
) -> BatchMapResponse:
    """
    Run a single MCP tool over a list of inputs in batch. Use this tool when you need to
    apply the same operation to multiple molecules efficiently in a single call.

    Note: MCP limits the size of parameter arrays to 128 items. For larger datasets, split into multiple batch_map calls.

    Example: To calculate molecular weight for 3 molecules, call batch_map with:
        tool_name="MolWt"
        inputs=[{"smiles": "CCO"}, {"smiles": "CC(=O)O"}, {"smiles": "C1=CC=CC=C1"}]

    Args:
        tool_name: Name of the MCP tool to call for each input (e.g., "MolWt", "NumRotatableBonds")
        inputs: List of argument dictionaries to map over (each dict is passed as tool arguments)
        ctx: MCP context (automatically injected)
        concurrency: Max concurrent calls (1-100). If None, auto-tunes based on input size.
        fail_fast: Stop on first error (default False)
        include_input: Include original input in each result item (default True)

    Returns:
        BatchMapResponse with results for each input
    """
    if concurrency is None:
        concurrency = optimal_concurrency(len(inputs))
    elif not 1 <= concurrency <= 100:
        raise ValueError("concurrency must be between 1 and 100")

    # Resolve tool name (strip prefix if needed, handle ambiguity)
    sample_input = inputs[0] if inputs else None
    actual_tool_name = await resolve_tool_name(tool_name, ctx, sample_input)

    sem = asyncio.Semaphore(concurrency)

    async def run_one(x: Dict[str, Any]) -> BatchItemResult:
        async with sem:
            try:
                out = await ctx.fastmcp.call_tool(actual_tool_name, x)
                return BatchItemResult(ok=True, input=x if include_input else None, output=out)
            except Exception as e:
                return BatchItemResult(ok=False, input=x if include_input else None, error=str(e))

    results: List[BatchItemResult] = []
    for coro in asyncio.as_completed([run_one(x) for x in inputs]):
        item = await coro
        results.append(item)
        if fail_fast and not item.ok:
            break

    return BatchMapResponse(tool_name=tool_name, concurrency=concurrency, results=results)
