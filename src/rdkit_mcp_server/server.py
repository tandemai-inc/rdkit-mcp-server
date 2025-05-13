import itertools
import os
import logging
from mcp.server.fastmcp import FastMCP

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

mcp = FastMCP("RDKit-MCP Server")

from . import tools as base_tools
from .rdkit.Chem.Descriptors import tools as rdkit_tools

# Modules to search for tools
TOOL_MODULES = [
    base_tools,
    rdkit_tools,
]


def get_tools_from_module(module):
    """Get all functions from the module that are decorated with @mcp.tool()"""
    return [
        getattr(module, func) for func in dir(module) if callable(getattr(module, func)) and not func.startswith("_")
    ]


def get_all_tools(tool_modules):
    tool_fns = itertools.chain.from_iterable(
        get_tools_from_module(module) for module in tool_modules
    )
    return tool_fns


def main():
    """Main function to run the MCP server."""
    # Determine transport method (default to stdio)
    transport = os.getenv("MCP_TRANSPORT", "sse").lower()
    host = os.getenv("MCP_HOST", "127.0.0.1") # Default host for SSE
    port_str = os.getenv("MCP_PORT", "8000")  # Default port for SSE

    try:
        port = int(port_str)
    except ValueError:
        logger.warning(f"Invalid MCP_PORT value '{port_str}'. Using default port 8000.")
        port = 8000
    
    logger.info("Registering tools with MCP server...")
    for tool in get_all_tools():
        try:
            mcp.add_tool(tool)
        except Exception as e:
            logger.error(f"Failed to register tool {tool.__name__}: {e}")

    logger.info(f"Starting RDKit MCP Server with transport: {transport}")
    if transport == "sse":
        logger.info(f"SSE transport selected. Listening on {host}:{port}")
        mcp.run(transport="sse")
    elif transport == "stdio":
        mcp.run(transport="stdio")
    else:
        logger.error(f"Unsupported transport type: {transport}. Defaulting to stdio.")
        mcp.run(transport="stdio")

# Note: The if __name__ == "__main__": block is typically placed in __main__.py
# for running the package directly with `python -m`.
# However, having a main() function here allows server.py to be run directly
# for testing if needed, or called from __main__.py.
