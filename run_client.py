import argparse
import asyncio

from rdkit_mcp_clients.openai import main as openai_client_main


parser = argparse.ArgumentParser(description="Run the OpenAI client with a prompt.")
parser.add_argument("--prompt", type=str, required=False, help="Prompt to pass to the client")
parser.add_argument("--model", "-m", type=str, default="gpt-4o-mini", help="OPENAI model")
args = parser.parse_args()

if __name__ == "__main__":
    asyncio.run(openai_client_main(args.prompt, args.model))
