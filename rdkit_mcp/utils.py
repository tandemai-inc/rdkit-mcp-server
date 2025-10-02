import base64
import logging
import pickle
from anyio import Path
from rdkit import Chem
from typing import Callable

from .types import EncodedFileModel, PickledMol, EncodedFile


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


def decode_mol(b64_pickled_mol: str) -> Chem.Mol:
    """
    Decode a base64 encoded pickled RDKit Mol object.

    Args:
        pickled_mol (str): Base64 encoded string of a pickled RDKit Mol object.

    Returns:
        Chem.Mol: The decoded RDKit Mol object.
    """
    pkl_data = base64.b64decode(b64_pickled_mol)
    mol = pickle.loads(pkl_data)
    if mol is None:
        raise ValueError("Failed to decode the molecule from the provided pickled data.")
    return mol


def encode_mol(mol: Chem.Mol) -> PickledMol:
    """
    Encode an RDKit Mol object into a base64 encoded pickled string.

    Args:
        mol (Chem.Mol): The RDKit Mol object to encode.

    Returns:
        str: Base64 encoded string of the pickled RDKit Mol object.
    """
    pkl_mol_data = pickle.dumps(mol)
    encoded_mol = base64.b64encode(pkl_mol_data).decode('utf-8')
    return encoded_mol


def encode_file_contents(file_path: Path, filename=None) -> EncodedFile:
    """
    Encode a file into a base64 encoded string.

    Args:
        file_path (Path): The path to the file to encode.
        filename (str, optional): A custom filename to use in the EncodedFileModel. Defaults to the name of the file at file_path.

    Returns:
        str: Base64 encoded string of the file contents.
    """
    with open(file_path, "rb") as f:
        file_data = f.read()

    if filename is None:
        filename = Path(file_path).name
    encoded_file_contents = base64.b64encode(file_data).decode('utf-8')
    encoded_file = EncodedFileModel(
        filename=filename,
        content=encoded_file_contents
    )
    return encoded_file


def decode_file_contents(b64_file_contents: str) -> str:
    """
    Decode a base64 encoded file.

    Args:
        b64_file_contents (str): The base64 encoded file contents.

    Returns:
        str: Contents of the file as a string.
    """
    file_data: bytes = base64.b64decode(b64_file_contents)
    return file_data.decode('utf-8')
