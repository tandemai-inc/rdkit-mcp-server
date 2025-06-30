import asyncio
import logging

from server.run import main

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


if __name__ == "__main__":
    # Configure logging for the script
    logger.info("Starting the RDKit MCP Server...")

    # Start the MCP server
    asyncio.run(main())
