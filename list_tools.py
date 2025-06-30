import argparse
import asyncio
import yaml

from rdkit_mcp.register_tools import register_tools
from app.settings import AppSettings, create_app_settings, get_app_settings
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


def load_settings(settings_path) -> AppSettings:
    if settings_path:
        with open(settings_path, "r") as f:
            yaml_data = yaml.safe_load(f)
            return create_app_settings(yaml_data)

    return get_app_settings()


async def list_tools():
    """Prints the module path of all tools registered to the MCP server."""
    args = parse_args()
    settings: AppSettings = load_settings(args.settings)
    print("Registered tools:")
    allow_list = settings.ALLOW_LIST
    block_list = settings.BLOCK_LIST
    await register_tools(mcp, allow_list=allow_list, block_list=block_list)
    tool_list = await mcp.list_tools()
    output = []
    for tool in tool_list:
        module_path = tool.annotations.module.title
        output.append(module_path)
    return output


if __name__ == "__main__":
    args = parse_args()
    settings = load_settings(args.settings)
    tool_list = asyncio.run(list_tools())
    print(f"Registered Tools: {len(tool_list)}")
    for tool in tool_list:
        print(f"- {tool}")
