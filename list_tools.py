import argparse
import asyncio
import yaml

from register_tools import register_tools
from run_server import mcp


def parse_args():
    parser = argparse.ArgumentParser(description="List registered tools from MCP server.")
    parser.add_argument(
        "--settings",
        type=str,
        required=False,
        help="Path to a YAML settings file."
    )
    return parser.parse_args()


def load_settings(settings_path):
    if settings_path:
        with open(settings_path, "r") as f:
            return yaml.safe_load(f)
    return {}


args = parse_args()
settings = load_settings(args.settings)


async def list_tools():
    """Prints the module path of all tools registered to the MCP server."""
    print("Registered tools:")
    allow_list = settings.get("allow_list", [])
    block_list = settings.get("block_list", [])
    await register_tools(mcp, allow_list=allow_list, block_list=block_list)
    tool_list = await mcp.list_tools()
    for tool in tool_list:
        module_path = tool.annotations.module.title
        print(f"- {module_path}")

if __name__ == "__main__":
    asyncio.run(list_tools())
