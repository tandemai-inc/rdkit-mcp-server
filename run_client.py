import argparse
import asyncio

from clients.openai import main as openai_client_main


parser = argparse.ArgumentParser(description="Run the OpenAI client with a prompt.")
parser.add_argument("--prompt", type=str, required=False, help="Prompt to pass to the client")
args = parser.parse_args()

if __name__ == "__main__":
    asyncio.run(openai_client_main(args.prompt))
