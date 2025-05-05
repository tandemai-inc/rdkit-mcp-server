import os
import logging
from mcp.server.fastmcp import FastMCP

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Initialize FastMCP server instance
# The server name might be used by clients to identify this server.
mcp = FastMCP("RDKit-MCP Server")

# Import tool functions from the tools module.
# The @mcp.tool() decorator in tools.py will register them with this 'mcp' instance.
logger.info("Importing RDKit tools...")
try:
    from . import tools
    logger.info("Successfully imported tools.")
except ImportError as e:
    logger.error(f"Failed to import tools module: {e}")
    # Depending on the desired behavior, we might exit or continue without tools.
    # For now, we'll log the error and continue, but the server might not have tools.
except Exception as e:
    logger.error(f"An unexpected error occurred during tool import: {e}")


def main():
    """Main function to run the MCP server."""
    # Determine transport method (default to stdio)
    transport = os.getenv("MCP_TRANSPORT", "stdio").lower()
    host = os.getenv("MCP_HOST", "127.0.0.1") # Default host for SSE
    port_str = os.getenv("MCP_PORT", "8000")  # Default port for SSE

    try:
        port = int(port_str)
    except ValueError:
        logger.warning(f"Invalid MCP_PORT value '{port_str}'. Using default port 8000.")
        port = 8000

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
