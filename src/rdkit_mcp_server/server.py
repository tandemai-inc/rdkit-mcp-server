import os
import logging
from mcp.server.fastmcp import FastMCP

from .tools.register_tools import register_tools

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


mcp = FastMCP("RDKit-MCP Server")


async def main():
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
    # TODO: Add settings to have whitelist and blacklist
    whitelist = []
    blacklist = []
    await register_tools(mcp, whitelist=whitelist, blacklist=blacklist)

    logger.info(f"Starting RDKit MCP Server with transport: {transport}")
    if transport == "sse":
        logger.info(f"SSE transport selected. Listening on {host}:{port}")
        await mcp.run_sse_async()
    elif transport == "stdio":
        mcp.run_stdio_async()
    else:
        logger.error(f"Unsupported transport type: {transport}. Defaulting to sse.")
        mcp.run_sse_async()
