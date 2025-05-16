import logging
from mcp.server.fastmcp import FastMCP
from typing import List, Callable

from . import base_tools
from .decorators import singleton
from .rdkit.Chem import Descriptors as rdkit_chem_descriptors

logger = logging.getLogger(__name__)


__all__ = ["register_tools", "get_tool_registry"]


# Modules to search for tools
TOOL_MODULES = [
    base_tools,
    rdkit_chem_descriptors,
]


def is_rdkit_tool(func):
    """Check if a function is decorated with @rdkit_tool."""
    # Access the original function in case of multiple wrappers
    while hasattr(func, "__wrapped__"):
        func = func.__wrapped__
    return hasattr(func, "_is_rdkit_tool")


def get_all_tools(tool_modules):
    for module in tool_modules:
        tool_iter = (
            getattr(module, func)
            for func in dir(module)
            # Only include if the function is decorated with rdkit_tool
            if is_rdkit_tool(getattr(module, func))
        )
        yield from tool_iter


def register_tools(mcp: FastMCP, whitelist: List[str]=None, blacklist: List[str]=None) -> None:
    """
    Register tools with the MCP server.

    Parameters:
    - mcp (FastMCP): The MCP server instance.
    """
    breakpoint()
    filter_list = ''
    if whitelist:
        filter_list = 'whitelist'
    elif blacklist:
        filter_list = 'blacklist'
    if whitelist and blacklist:
        logger.warning("Both whitelist and blacklist provided. Using whitelist.")

    for tool in get_all_tools(TOOL_MODULES):
        try:
            tool_name = f'{tool.__module__}.{tool.__name__}'
            if filter_list == 'whitelist' and tool_name not in whitelist:
                logger.warning(f"Tool {tool_name} not in whitelist. Skipping.")
                continue
            if filter_list == 'blacklist' and tool_name in blacklist:
                logger.warning(f"Tool {tool_name} in blacklist. Skipping.")
                continue
            mcp.add_tool(
                tool,
                name=tool_name,	
                description=tool.__doc__,
            )
        except Exception as e:
            logger.error(f"Failed to register tool {tool.__name__}: {e}")

@singleton
class ToolRegistry:

    def __init__(self):
        self._registry = {}

    def register_tool(self, name: str, func: Callable, enabled: bool = True):
        """Register a tool and set its enabled state."""
        self._registry[name] = {
            "func": func,
            "enabled": enabled
        }

    def set_tool_enabled(self, name: str, enabled: bool):
        """Enable or disable a tool by name."""
        if name in self._registry:
            self._registry[name]["enabled"] = enabled

    def is_tool_enabled(self, name: str) -> bool:
        """Check if a tool is enabled."""
        return self._registry.get(name, {}).get("enabled", False)

    def get_tool(self, name: str) -> Callable:
        """Get the tool function, enforcing the enabled state."""
        tool_info = self._registry.get(name)
        if tool_info and tool_info["enabled"]:
            return tool_info["func"]
        raise RuntimeError(f"Tool '{name}' is disabled or not found.")

    def list_tools(self):
        """List all registered tools and their states."""
        return {name: info["enabled"] for name, info in self._registry.items()}


# Singleton instance accessor
def get_tool_registry() -> ToolRegistry:
    return ToolRegistry()