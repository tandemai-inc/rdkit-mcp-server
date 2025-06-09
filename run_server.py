import asyncio
import logging
from rdkit_mcp.server import main, logger

if __name__ == "__main__":
    # Configure logging for the script
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger.info("Starting the RDKit MCP Server...")

    # Start the MCP server
    asyncio.run(main())
