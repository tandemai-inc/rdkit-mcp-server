import itertools
import os
import logging
from mcp.server.fastmcp import FastMCP

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

mcp = FastMCP("RDKit-MCP Server")

from . import tools as base_tools
from .rdkit.Chem import Descriptors as descriptor_tools

# Modules to search for tools
TOOL_MODULES = [
    base_tools,
    descriptor_tools,
]

def get_all_tools(tool_modules):
    for module in tool_modules:
        tool_iter = (
            getattr(module, func)
            for func in dir(module)
            if callable(getattr(module, func))
            # Filter out private functions and those not starting with "_"
            and not func.startswith("_")
            # Only include functions defined in the module
            and getattr(getattr(module, func), "__module__", None) == module.__name__
        )
        yield from tool_iter


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
    for tool in get_all_tools(TOOL_MODULES):
        try:
            mcp.add_tool(
                tool,
                name=tool.__name__,
                description=tool.__doc__,
            )
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
