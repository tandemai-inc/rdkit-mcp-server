from typing import Callable
import logging


logger = logging.getLogger(__name__)


def is_rdkit_tool(func: Callable) -> bool:
    """Check if a function is decorated with @rdkit_tool."""
    # Access the original function in case of multiple wrappers
    while hasattr(func, "__wrapped__"):
        func = func.__wrapped__
    return hasattr(func, "_is_rdkit_tool")


def singleton(cls):
    """Add as a decorator to a class to make it a singleton."""
    instances = {}

    def get_instance(*args, **kwargs):
        if cls not in instances:
            logger.debug("Creating singleton instance of %s with args: %s, kwargs: %s", cls.__name__, args, kwargs)
            instances[cls] = cls(*args, **kwargs)
        elif kwargs and instances[cls] is not None:
            # If kwargs are provided, create a new instance
            # This fixes issue where ToolSettings can be instantiated before yaml is loaded,
            # and then we can never get the correct settings.
            logger.debug("Overwriting singleton instance of %s with kwargs: %s", cls.__name__, kwargs)
            instances[cls] = cls(*args, **kwargs)
        return instances[cls]
    return get_instance
