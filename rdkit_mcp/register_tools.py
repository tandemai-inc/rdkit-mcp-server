import logging
from mcp.server.fastmcp import FastMCP
from typing import Callable, Iterable, List
from rdkit_mcp import get_rdkit_tools

logger = logging.getLogger(__name__)


__all__ = ["register_tools"]


async def register_tools(mcp: FastMCP, allow_list: List[str] = None, block_list: List[str] = None) -> None:
    """
    Register tools with the MCP server.

    Parameters:
    - mcp (FastMCP): The MCP server instance.
    """
    filter_list = ''
    if allow_list:
        filter_list = 'allow_list'
    elif block_list:
        filter_list = 'block_list'
    if allow_list and block_list:
        logger.warning("Both allow_list and block_list of tools provided. Using allow_list.")

    rdkit_tool_iter: Iterable[Callable] = get_rdkit_tools()
    # Loop through all tools and register them with the MCP server
    for tool_fn in rdkit_tool_iter:
        try:
            tool_module = tool_fn.tool_annotations['module'].title
            tool_name = tool_fn.tool_name or tool_fn.__name__
            if not tool_fn.tool_enabled:
                logger.debug(f"Tool {tool_name} is disabled. Skipping.")
                continue
            if filter_list == 'allow_list' and not _tool_module_matches(tool_module, allow_list):
                logger.debug(f"Tool {tool_name} not matched by allow_list. Skipping.")
                continue
            if filter_list == 'block_list' and _tool_module_matches(tool_module, block_list):
                logger.debug(f"Tool {tool_name} matched by block_list. Skipping.")
                continue
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


def _tool_module_matches(name: str, patterns: List[str]) -> bool:
    """
    Returns True if any pattern in patterns is a substring of name.
    """
    if not patterns:
        return False
    return any(pattern in name for pattern in patterns)
