import logging
from mcp.server.fastmcp import FastMCP
from typing import Callable, Iterable, List
from rdkit_mcp import get_rdkit_tools

logger = logging.getLogger(__name__)


__all__ = ["register_tools", "collect_tools"]


def _tool_module_matches(tool_module: str, filter_list: List[str]) -> bool:
    """
    Returns True if any pattern in filter_list is a substring of name.
    """
    tool_module_fn_name = tool_module.split('.')[-1]

    match_found = False
    for pattern in filter_list:
        # Check for straight match of the function name.
        if tool_module in filter_list:
            match_found = True
            break
        # Check if the tool_module is an exact match with the last part of the module path.
        # Example: if tool_module is 'MolWt' it should match rdkit_mcp.Chem.Descriptors.MolWt,
        # but not `rdkit_mcp.Chem.Descriptors.ExactMolWt`
        pattern_fn_name = pattern.split('.')[-1]
        if tool_module_fn_name == pattern_fn_name:
            match_found = True
            break
        # Check if the tool is in a module that is being filtered
        # For Example, `rdkit_mcp.Chem.rdMolDescriptors` in the filter_list
        # should match `rdkit_mcp.Chem.rdMolDescriptors.CalcPBF`
        if tool_module.startswith(pattern):
            match_found = True
            break
    return match_found


def collect_tools(allow_list: List[str] = None, block_list: List[str] = None) -> Iterable[Callable]:
    """
    Collect tools based on allow_list and block_list.

    Parameters:
    - allow_list (List[str]): List of tool names to allow.
    - block_list (List[str]): List of tool names to block.

    Returns:
    - Iterable[Callable]: An iterable of tool functions.
    """
    rdkit_tool_iter: Iterable[Callable] = get_rdkit_tools()
    filter_list = None
    if allow_list:
        filter_list = 'allow_list'
    elif block_list:
        filter_list = 'block_list'
    if allow_list and block_list:
        logger.warning("Both allow_list and block_list of tools provided. Using allow_list.")

    rdkit_tool_iter: Iterable[Callable] = get_rdkit_tools()
    for tool_fn in rdkit_tool_iter:
        tool_module = tool_fn.tool_annotations['module'].title
        if not tool_fn.tool_enabled:
            logger.debug(f"Tool {tool_fn.__name__} is disabled. Skipping.")
            continue
        if filter_list == 'allow_list' and not _tool_module_matches(tool_module, allow_list):
            logger.debug(f"Tool {tool_fn.__name__} not matched by allow_list. Skipping.")
            continue
        if filter_list == 'block_list' and _tool_module_matches(tool_module, block_list):
            logger.debug(f"Tool {tool_fn.__name__} matched by block_list. Skipping.")
            continue
        yield tool_fn


async def register_tools(mcp: FastMCP, allow_list: List[str] = None, block_list: List[str] = None) -> None:
    """
    Register tools with the MCP server.

    Parameters:
    - mcp (FastMCP): The MCP server instance.
    """
    # Loop through all tools and register them with the MCP server
    filtered_tool_list = collect_tools(allow_list=allow_list, block_list=block_list)
    for tool_fn in filtered_tool_list:
        try:
            tool_name = tool_fn.tool_name or tool_fn.__name__
            # Add tool the MCP Server
            # These properties on the function are set by the rdkit_tool decorator
            tool_description = getattr(tool_fn, 'tool_description', tool_fn.__doc__)
            tool_annotations = getattr(tool_fn, 'tool_annotations', None)
            mcp.add_tool(
                tool_fn,
                name=tool_name,
                description=tool_description,
                annotations=tool_annotations,
            )
        except Exception as e:
            logger.error(f"Failed to register tool {tool_fn.__name__}: {e}")
    tool_count = len(await mcp.list_tools())
    if tool_count == 0:
        raise RuntimeError("No tools registered with MCP server. Please check your allow_lists/block_lists.")
    logger.info(f"Registered {tool_count} tools with MCP server.")
