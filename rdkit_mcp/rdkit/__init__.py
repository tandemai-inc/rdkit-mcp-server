from ..tools.utils import is_rdkit_tool
from typing import Iterable, Callable
import importlib
import pkgutil
import os


def get_rdkit_tools() -> Iterable[Callable]:
    """Walk packages in the rdkit_mcp.rdkit package and yield all callable rdkit tools."""
    pkg_dir = os.path.dirname(__file__)
    pkg_name = __package__ or "rdkit_mcp.rdkit"

    # Recursively walk through all modules in the current package
    for finder, name, ispkg in pkgutil.walk_packages([pkg_dir], prefix=pkg_name + "."):
        try:
            module = importlib.import_module(name)
        except Exception:
            continue  # Skip modules that fail to import

        for attr_name in dir(module):
            attr = getattr(module, attr_name)
            if is_rdkit_tool(attr) and getattr(attr, "tool_enabled", True):
                yield attr
