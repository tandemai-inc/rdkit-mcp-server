import functools
from functools import wraps
from typing import Type, TypeVar

T = TypeVar('T')


def singleton(cls: Type[T]) -> Type[T]:
    """A decorator to make a class a Singleton."""
    instances = {}

    @wraps(cls)
    def get_instance(*args, **kwargs) -> T:
        if cls not in instances:
            instances[cls] = cls(*args, **kwargs)
        return instances[cls]
    return get_instance


def rdkit_tool(func):
    """Decorator to mark a function as an RDKit tool."""
    # Mark the function with a custom attribute
    func._is_rdkit_tool = True
    return func
