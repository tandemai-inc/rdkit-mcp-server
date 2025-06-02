import argparse
import logging
from mcp.server.fastmcp import FastMCP

from .tools.register_tools import register_tools

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="RDKit MCP Server")
parser.add_argument("--port", type=int, help="Port to run the server on", default=8000)
parser.add_argument("--transport", choices=["sse", "stdio"], help="Transport method (sse or stdio)", default="sse")
parser.add_argument("--host", type=str, help="Host to run the server on", default="127.0.0.1")

mcp = FastMCP("RDKit-MCP Server")


async def main():
    """Main function to run the MCP server."""
    args, _ = parser.parse_known_args()
    transport = args.transport

    # TODO: Add settings to have whitelist and blacklist
    whitelist = []
    blacklist = []
    logger.info("Registering tools with MCP server...")
    await register_tools(mcp, whitelist=whitelist, blacklist=blacklist)

    logger.info(f"Starting RDKit MCP Server with transport: {transport}")
    if transport == "sse":
        logger.info(f"Server running on {args.host}:{args.port} using SSE transport.")
        mcp.settings.host = args.host
        mcp.settings.port = args.port
        await mcp.run_sse_async()
    elif transport == "stdio":
        logger.info("Server running using stdio transport.")
        await mcp.run_stdio_async()
