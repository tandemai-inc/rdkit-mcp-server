import asyncio
import argparse
import logging
import os
import yaml

from fastapi import Request, Response
from fastapi.responses import FileResponse, JSONResponse
from mcp.server.fastmcp import FastMCP

from rdkit_mcp.register_tools import register_tools
from rdkit_mcp.settings import ToolSettings

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(name)s - %(message)s'
)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="RDKit MCP Server")
parser.add_argument("--port", type=int, help="Port to run the server on", default=8000)
parser.add_argument("--transport", choices=["sse", "stdio"], help="Transport method (sse or stdio)", default="sse")
parser.add_argument("--host", type=str, help="Host to run the server on", default="127.0.0.1")
parser.add_argument("--settings", type=str, help="Path to YAML settings file", default=None)
parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (debug) logging")

mcp = FastMCP("RDKit-MCP Server")


@mcp.custom_route("/files", methods=["GET"])
async def list_files(request: Request):
    settings = ToolSettings()
    file_dir = str(settings.FILE_DIR)
    try:
        files = os.listdir(file_dir)
        return JSONResponse({"files": files})
    except Exception as e:
        return JSONResponse({"error": str(e)})


@mcp.custom_route("/files/{filename}", methods=["GET"])
async def get_file(request: Request):
    filename = request.path_params.get("filename")
    settings = ToolSettings()
    file_dir = str(settings.FILE_DIR)
    file_path = os.path.join(file_dir, filename)
    if not os.path.isfile(file_path):
        return JSONResponse(status_code=404, content={"error": "File not found"})
    return FileResponse(path=file_path, filename=filename)    


async def main():
    """Main function to run the MCP server."""
    args, _ = parser.parse_known_args()
    transport = args.transport

    # Load settings from YAML file if provided
    settings_data = {}
    if args.settings:
        try:
            with open(args.settings, "r") as f:
                settings_data = yaml.safe_load(f)
            logger.info(f"Loaded settings from {args.settings}: {settings_data}")
        except FileNotFoundError:
            logger.error(f"Settings file not found: {args.settings}")
        except yaml.YAMLError as e:
            logger.error(f"Error parsing YAML settings file: {e}")
    settings: ToolSettings = ToolSettings(**settings_data)

    # Register tools with the MCP server
    allow_list = settings.ALLOW_LIST
    block_list = settings.BLOCK_LIST
    logger.info("Registering tools with MCP server...")
    await register_tools(mcp, allow_list=allow_list, block_list=block_list)

    # Start server
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
    args, _ = parser.parse_known_args()
    if args.verbose:
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled.")

    logger.info("Starting the RDKit MCP Server...")

    # Start the MCP server
    asyncio.run(main())
