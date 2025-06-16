import logging
import os
from typing import Union
from mcp.server.fastmcp.exceptions import ToolError
from pathlib import Path
from rdkit import Chem

from .utils import rdkit_tool, OUTPUT_DIR
from .types import Smiles

logger = logging.getLogger(__name__)


@rdkit_tool(enabled=False)
def smiles_to_sdf(smiles: Smiles) -> Path:
    """
    Converts a SMILES string to an SDF file.

    Args:
        smiles: The SMILES representation of the molecule.
    Returns:
        An SDF string representation of the molecule.
    """
    logger.info(f"Converting SMILES to SDF: {smiles[:30]}...")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ToolError(f"Invalid or unparsable SMILES string: {smiles}")

    sdf_string = Chem.MolToMolBlock(mol)
    # Write to SDF file
    filename = f"{Chem.MolToSmiles(mol)}.sdf"
    output_path = Path(os.path.join(OUTPUT_DIR, filename))
    with open(output_path, "w") as f:
        f.write(sdf_string)
    return output_path


@rdkit_tool(enabled=False)
def sdf_to_smiles(sdf_path: Union[str, Path]) -> Smiles:
    """
    Converts an SDF file to a SMILES string.

    Args:
        sdf_path: The path to the SDF file.
    Returns:
        The SMILES representation of the molecule.
    """
    logger.info(f"Converting SDF to SMILES: {sdf_path}")
    if isinstance(sdf_path, str):
        sdf_path = Path(sdf_path)

    if not sdf_path.exists():
        raise ToolError(f"SDF file does not exist: {sdf_path}")

    suppl = Chem.SDMolSupplier(str(sdf_path))
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ToolError(f"Failed to read any valid molecule from SDF file: {sdf_path}")

    smiles = Chem.MolToSmiles(mol)
    return smiles


def get_base_tools():
    """
    Get all base tools defined in this module.

    Returns:
        A list of tool functions.
    """
    return [
        smiles_to_sdf,
        sdf_to_smiles,
    ]
