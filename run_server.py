import asyncio
import argparse
import logging
import yaml
from mcp.server.fastmcp import FastMCP

from rdkit_mcp.register_tools import register_tools
from settings import AppSettings, create_app_settings, get_app_settings

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="RDKit MCP Server")
parser.add_argument("--port", type=int, help="Port to run the server on", default=8000)
parser.add_argument("--transport", choices=["sse", "stdio"], help="Transport method (sse or stdio)", default="sse")
parser.add_argument("--host", type=str, help="Host to run the server on", default="127.0.0.1")
parser.add_argument("--settings", type=str, help="Path to YAML settings file", default=None)

mcp = FastMCP("RDKit-MCP Server")


async def main():
    """Main function to run the MCP server."""
    args, _ = parser.parse_known_args()
    transport = args.transport
    if not args.settings:
        logger.debug("No settings file provided, using default settings.")
        settings: AppSettings = get_app_settings()
    else:
        with open(args.settings, "r") as f:
            yaml_settings = yaml.safe_load(f)
        settings: AppSettings = create_app_settings(yaml_settings)
        logger.info(f"Loaded settings from {args.settings}: {yaml_settings}")

    allow_list = settings.ALLOW_LIST
    block_list = settings.BLOCK_LIST
    logger.info("Registering tools with MCP server...")
    await register_tools(mcp, allow_list=allow_list, block_list=block_list)
    logger.info(f"Starting RDKit MCP Server with transport: {transport}")
    if transport == "sse":
        logger.info(f"Server running on {args.host}:{args.port} using SSE transport.")
        mcp.settings.host = args.host
        mcp.settings.port = args.port
        await mcp.run_sse_async()
    elif transport == "stdio":
        logger.info("Server running using stdio transport.")
        await mcp.run_stdio_async()


if __name__ == "__main__":
    # Configure logging for the script
    logger.info("Starting the RDKit MCP Server...")

    # Start the MCP server
    asyncio.run(main())
