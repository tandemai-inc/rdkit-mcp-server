from .Chem import Descriptors, AllChem, rdMMPA
from ..utils import is_rdkit_tool
from typing import Iterable, Callable


# Add new modules wrapped in the MCP server here
TOOL_MODULES = [
    Descriptors,
    AllChem,
    rdMMPA,
]


def get_rdkit_tools() -> Iterable[Callable]:
    for module in TOOL_MODULES:
        tool_iter = (
            getattr(module, func)
            for func in dir(module)
            if is_rdkit_tool(getattr(module, func))
            # Tool is not disabled
            and not getattr(getattr(module, func), "tool_disabled", False)
        )
        yield from tool_iter
