import sys
from .server import main, logger

if __name__ == "__main__":
    logger.info(f"Running RDKit MCP Server via __main__.py (args: {sys.argv})")
    main()
