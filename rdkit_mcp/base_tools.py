import logging
import os
from typing import Union
from mcp.server.fastmcp.exceptions import ToolError
from pathlib import Path
from rdkit import Chem

from rdkit_mcp.settings import ToolSettings
from .decorators import rdkit_tool
from .types import PickledMol, Smiles, EncodedFile
from .utils import encode_mol, decode_mol, encode_file_contents, decode_file_contents

logger = logging.getLogger(__name__)


@rdkit_tool()
def smiles_to_mol(smiles: Smiles) -> PickledMol:
    """
    Converts a SMILES string into a base 64 encoded pickled RDKit Mol object.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ToolError(f"Invalid or unparsable SMILES string: {smiles}")
    encoded_mol = encode_mol(mol)
    return encoded_mol


@rdkit_tool(description="Converts a pickled RDKit mol object to a SMILES string.")
def mol_to_smiles(pmol: PickledMol) -> Smiles:
    """
    Converts a pickled RDKit mol object to a SMILES string.
    """
    mol = decode_mol(pmol)
    if mol is None:
        raise ToolError(f"Failed to decode the pickled RDKit Mol object.")
    smiles = Chem.MolToSmiles(mol)
    return smiles


@rdkit_tool(enabled=False)
def smiles_to_sdf(smiles: Smiles) -> EncodedFile:
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

    return encode_file_contents(sdf_string)


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


@rdkit_tool(description="Writes a pickled RDKit Mol object to an SDF file and returns the file path.")
def mol_to_sdf(pmol: PickledMol, filename=None) -> EncodedFile:
    """
    Writes a pickled RDKit Mol object to an SDF file.

    Args:
        pmol: The pickled and base64-encoded RDKit Mol object.
    Returns:
        The path to the written SDF file.
    """
    mol = decode_mol(pmol)
    if mol is None:
        raise ToolError("Failed to decode the pickled RDKit Mol object.")

    sdf_string = Chem.MolToMolBlock(mol)

    if filename is None:
        filename = f"{Chem.MolToSmiles(mol)}.sdf"
    if not filename.endswith('.sdf'):
        filename += '.sdf'
    return encode_file_contents(sdf_string)


@rdkit_tool(description="Converts a PDB file to a pickled RDKit Mol object.")
def pdb_to_mol(pdb_path: Union[str, Path]) -> PickledMol:
    """
    Converts a PDB file to a pickled RDKit Mol object.

    Args:
        pdb_path: The path to the PDB file.
    Returns:
        A base64 encoded pickled RDKit Mol object.
    """
    if isinstance(pdb_path, str):
        pdb_path = Path(pdb_path)

    if not pdb_path.exists():
        raise ToolError(f"PDB file does not exist: {pdb_path}")

    mol = Chem.MolFromPDBFile(str(pdb_path))
    if mol is None:
        raise ToolError(f"Failed to read molecule from PDB file: {pdb_path}")

    encoded_mol = encode_mol(mol)
    return encoded_mol


@rdkit_tool(description="Converts a pickled RDKit Mol object to a PDB file and returns the file path.")
def mol_to_pdb(pmol: PickledMol, filename=None) -> EncodedFile:
    """
    Converts a pickled RDKit Mol object to a PDB file.

    Args:
        pmol: The pickled and base64-encoded RDKit Mol object.
    Returns:
        The path to the written PDB file.
    """
    mol = decode_mol(pmol)
    if mol is None:
        raise ToolError("Failed to decode the pickled RDKit Mol object.")

    pdb_string = Chem.MolToPDBBlock(mol)

    if filename is None:
        filename = f"{Chem.MolToSmiles(mol)}.pdb"
    if not filename.endswith('.pdb'):
        filename += '.pdb'
    return encode_file_contents(pdb_string)


@rdkit_tool(description="Converts an SDF file to a pickled RDKit Mol object.")
def sdf_to_mol(sdf_path: Union[str, Path]) -> PickledMol:
    """
    Converts an SDF file to a pickled RDKit Mol object.

    Args:
        sdf_path: The path to the SDF file.
    Returns:
        A base64 encoded pickled RDKit Mol object.
    """
    if isinstance(sdf_path, str):
        sdf_path = Path(sdf_path)

    if not sdf_path.exists():
        raise ToolError(f"SDF file does not exist: {sdf_path}")

    suppl = Chem.SDMolSupplier(str(sdf_path))
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ToolError(f"Failed to read any valid molecule from SDF file: {sdf_path}")

    encoded_mol = encode_mol(mol)
    return encoded_mol


@rdkit_tool(description="Converts a pickled RDKit Mol object to an SDF file and returns the file path.")
def mol_to_sdf(pmol: PickledMol, filename=None) -> EncodedFile:
    """
    Converts a pickled RDKit Mol object to an SDF file.

    Args:
        pmol: The pickled and base64-encoded RDKit Mol object.
    Returns:
        The path to the written SDF file.
    """
    mol = decode_mol(pmol)
    if mol is None:
        raise ToolError("Failed to decode the pickled RDKit Mol object.")

    sdf_string = Chem.MolToMolBlock(mol)

    if filename is None:
        filename = f"{Chem.MolToSmiles(mol)}.sdf"
    if not filename.endswith('.sdf'):
        filename += '.sdf'
    return encode_file_contents(sdf_string)
