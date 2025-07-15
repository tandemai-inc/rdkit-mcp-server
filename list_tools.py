import argparse
import asyncio
import yaml
from typing import Callable, Iterable

from rdkit_mcp.register_tools import collect_tools
from rdkit_mcp.settings import ToolSettings


def parse_args():
    parser = argparse.ArgumentParser(description="List registered tools from MCP server.")
    parser.add_argument(
        "--settings",
        type=str,
        required=False,
        help="Path to a YAML settings file."
    )
    return parser.parse_args()


def load_settings(settings_path) -> ToolSettings:
    """Load settings from a YAML file or return default settings."""
    yaml_data = {}
    if settings_path:
        with open(settings_path, "r") as f:
            yaml_data = yaml.safe_load(f)
    return ToolSettings(**yaml_data)


async def list_tools(allow_list=None, block_list=None):
    """Returns the list of modules for all tools being registered to the MCP server."""
    allow_list = allow_list or []
    block_list = block_list or []
    tool_list: Iterable[Callable] = collect_tools(allow_list=allow_list, block_list=block_list)
    output = []
    for tool in tool_list:
        module_path = tool.annotations.module.title
        output.append(module_path)
    return output


if __name__ == "__main__":
    args = parse_args()
    settings: ToolSettings = load_settings(args.settings)
    allow_list = settings.ALLOW_LIST
    block_list = settings.BLOCK_LIST
    tool_list = asyncio.run(list_tools(allow_list=allow_list, block_list=block_list))
    print(f"Registered Tools: {len(tool_list)}")
    for tool in tool_list:
        print(f"- {tool}")
