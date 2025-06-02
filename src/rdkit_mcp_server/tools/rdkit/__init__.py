from .Chem import Descriptors, AllChem, rdMMPA, rdMolDescriptors
from ..utils import is_rdkit_tool
from typing import Iterable, Callable


# Add new modules wrapped in the MCP server here
TOOL_MODULES = [
    Descriptors,
    AllChem,
    rdMMPA,
    rdMolDescriptors,
]


def get_rdkit_tools() -> Iterable[Callable]:
    for module in TOOL_MODULES:
        tool_iter = (
            getattr(module, func)
            for func in dir(module)
            if is_rdkit_tool(getattr(module, func))
            # Tool is enabled
            and getattr(getattr(module, func), "tool_enabled", True)
        )
        yield from tool_iter
