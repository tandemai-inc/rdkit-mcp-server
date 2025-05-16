from  .base_tools import register_tools as base_tools
from .rdkit.Chem.Descriptors import register_tools as rdkit_chem_descriptor_tools
from mcp.server.fastmcp import FastMCP


__all__ = ["register_tools"]

# Modules to search for tools
tool_registry_fns = [
    base_tools,
    rdkit_chem_descriptor_tools,
]

def register_tools(mcp: FastMCP) -> None:
    """
    Register tools with the MCP server.

    Parameters:
    - mcp (FastMCP): The MCP server instance.
    """
    for tool_register_fn in tool_registry_fns:
        # Register tools from the specified modules
        tool_register_fn(mcp)
