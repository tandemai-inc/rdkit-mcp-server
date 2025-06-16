import logging
from mcp.server.fastmcp import FastMCP
from typing import Callable, Iterable, List
from rdkit_mcp import get_rdkit_tools

logger = logging.getLogger(__name__)


__all__ = ["register_tools"]


async def register_tools(mcp: FastMCP, whitelist: List[str] = None, blacklist: List[str] = None) -> None:
    """
    Register tools with the MCP server.

    Parameters:
    - mcp (FastMCP): The MCP server instance.
    """
    filter_list = ''
    if whitelist:
        filter_list = 'whitelist'
    elif blacklist:
        filter_list = 'blacklist'
    if whitelist and blacklist:
        logger.warning("Both whitelist and blacklist of tools provided. Using whitelist.")

    rdkit_tool_iter: Iterable[Callable] = get_rdkit_tools()

    # Loop through all tools and register them with the MCP server
    for tool_fn in rdkit_tool_iter:
        try:
            tool_name = getattr(tool_fn, 'tool_name', tool_fn.__name__)
            if not tool_fn.tool_enabled:
                logger.debug(f"Tool {tool_name} is disabled. Skipping.")
                continue
            if filter_list == 'whitelist' and tool_name not in whitelist:
                logger.debug(f"Tool {tool_name} not in whitelist. Skipping.")
                continue
            if filter_list == 'blacklist' and tool_name in blacklist:
                logger.debug(f"Tool {tool_name} in blacklist. Skipping.")
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
        raise RuntimeError("No tools registered with MCP server. Please check your whitelists/blacklists.")
    logger.info(f"Registered {tool_count} tools with MCP server.")
