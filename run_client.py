import asyncio
from src.clients.openai import main as openai_client_main


if __name__=="__main__":
    asyncio.run(openai_client_main())
