from typing import Callable


def rdkit_tool(func: Callable) -> Callable:
    """Decorator to mark a function as an RDKit tool."""
    # Mark the function with a custom attribute
    func._is_rdkit_tool = True
    return func


def is_rdkit_tool(func: Callable) -> bool:
    """Check if a function is decorated with @rdkit_tool."""
    # Access the original function in case of multiple wrappers
    while hasattr(func, "__wrapped__"):
        func = func.__wrapped__
    return hasattr(func, "_is_rdkit_tool")
