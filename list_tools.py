import argparse
import asyncio
import yaml

from rdkit_mcp.register_tools import collect_tools
from rdkit_mcp.settings import AppSettings, create_app_settings, get_app_settings
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


async def list_tools(allow_list=None, block_list=None):
    """Prints the module path of all tools registered to the MCP server."""
    allow_list = allow_list or []
    block_list = block_list or []
    tool_list = collect_tools(allow_list=allow_list, block_list=block_list)
    output = []
    for tool in tool_list:
        module_path = f'{tool.__module__}.{tool.__name__}'
        output.append(module_path)
    return output


if __name__ == "__main__":
    args = parse_args()
    settings: AppSettings = load_settings(args.settings)
    allow_list = settings.ALLOW_LIST
    block_list = settings.BLOCK_LIST
    tool_list = asyncio.run(list_tools(allow_list=allow_list, block_list=block_list))
    print(f"Registered Tools: {len(tool_list)}")
    for tool in tool_list:
        print(f"- {tool}")
