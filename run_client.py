import asyncio
from src.clients.mcp_client import main as mcp_main
from src.clients.openai import main as openai_main


if __name__=="__main__":
    OPENAI = 'openai'
    MCP = 'mcp'

    client_to_use = OPENAI

    if client_to_use == OPENAI:
        main_fn = openai_main
    else:
        main_fn = mcp_main

    asyncio.run(main_fn())
