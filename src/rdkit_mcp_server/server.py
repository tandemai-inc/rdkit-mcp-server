import os
import logging
from mcp.server.fastmcp import FastMCP

from .tools import register_tools

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


mcp = FastMCP("RDKit-MCP Server")


def main():
    """Main function to run the MCP server."""
    # Determine transport method (default to stdio)
    transport = os.getenv("MCP_TRANSPORT", "sse").lower()
    host = os.getenv("MCP_HOST", "127.0.0.1")  # Default host for SSE
    port_str = os.getenv("MCP_PORT", "8000")  # Default port for SSE

    try:
        port = int(port_str)
    except ValueError:
        logger.warning(f"Invalid MCP_PORT value '{port_str}'. Using default port 8000.")
        port = 8000

    logger.info("Registering tools with MCP server...")
    register_tools(mcp)

    logger.info(f"Starting RDKit MCP Server with transport: {transport}")
    if transport == "sse":
        logger.info(f"SSE transport selected. Listening on {host}:{port}")
        mcp.run(transport="sse")
    elif transport == "stdio":
        mcp.run(transport="stdio")
    else:
        logger.error(f"Unsupported transport type: {transport}. Defaulting to stdio.")
        mcp.run(transport="stdio")
