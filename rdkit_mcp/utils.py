import os
from typing import Callable

OUTPUT_DIR = os.path.join(os.getcwd(), 'outputs')


def is_rdkit_tool(func: Callable) -> bool:
    """Check if a function is decorated with @rdkit_tool."""
    # Access the original function in case of multiple wrappers
    while hasattr(func, "__wrapped__"):
        func = func.__wrapped__
    return hasattr(func, "_is_rdkit_tool")
