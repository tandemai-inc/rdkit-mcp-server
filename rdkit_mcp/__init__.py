from rdkit_mcp.utils import is_rdkit_tool
from typing import Iterable, Callable
import importlib
import pkgutil
import os
import logging

logger = logging.getLogger(__name__)


def get_rdkit_tools() -> Iterable[Callable]:
    """Walk packages in the rdkit_mcp.rdkit package and yield all callable rdkit tools."""
    pkg_dir = os.path.dirname(__file__)
    pkg_name = __package__ or "rdkit_mcp"
    # Recursively walk through all modules in the current package
    for _, name, _ in pkgutil.walk_packages(path=[pkg_dir], prefix=f"{pkg_name}."):
        try:
            module = importlib.import_module(name)
        except Exception as e:
            logger.debug(f"Failed to import module {name}: {e}")
            continue  # Skip modules that fail to import

        for attr_name in dir(module):
            attr = getattr(module, attr_name)
            if is_rdkit_tool(attr) and getattr(attr, "tool_enabled", True):
                yield attr
