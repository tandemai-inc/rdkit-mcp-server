from io import BytesIO
import logging
import tempfile
from typing import Union
from mcp.server.fastmcp.exceptions import ToolError
from pathlib import Path
from rdkit import Chem

from .decorators import rdkit_tool
from .types import PickledMol, Smiles, EncodedFile
from .utils import encode_mol, decode_mol, encode_file_contents

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


@rdkit_tool()
def mol_to_sdf(pmol: PickledMol, filename: Union[str, None] = None) -> EncodedFile:
    """
    Writes a pickled RDKit Mol object to an SDF file.

    Args:
        pmol: The pickled and base64-encoded RDKit Mol object.
        filename: Optional filename for the SDF file. If not provided, a name will be generated based on the SMILES string.
    Returns:
        A base64 encoded contents of an SDF file.
    """
    mol = decode_mol(pmol)
    if mol is None:
        raise ToolError("Failed to decode the pickled RDKit Mol object.")

    sdf_string = Chem.MolToMolBlock(mol)

    if filename is None:
        filename = f"{Chem.MolToSmiles(mol)}.sdf"
    if not filename.endswith('.sdf'):
        filename += '.sdf'

    with tempfile.NamedTemporaryFile(suffix=".sdf") as temp_sdf_file:
        temp_sdf_file.write(sdf_string.encode('utf-8'))
        temp_sdf_file.flush()
        return encode_file_contents(temp_sdf_file.name, filename=filename)


@rdkit_tool()
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


@rdkit_tool()
def pdb_contents_to_mol(pdb_contents: str) -> PickledMol:
    """
    Converts the contents of a PDB file to a pickled RDKit Mol object.

    Args:
        pdb_contents: The contents of the PDB file.
    Returns:
        A base64 encoded pickled RDKit Mol object.
    """
    mol = Chem.MolFromPDBBlock(pdb_contents)
    if mol is None:
        raise ToolError(f"Failed to read molecule from PDB contents.")
    encoded_mol = encode_mol(mol)
    return encoded_mol


@rdkit_tool()
def mol_to_pdb(pmol: PickledMol, filename: Union[str, None] = None) -> EncodedFile:
    """
    Converts a pickled RDKit Mol object to a PDB file.

    Args:
        pmol: The pickled and base64-encoded RDKit Mol object.
        filename: Optional filename for the PDB file.
    Returns:
       A base64 encoded contents of an PDB file.
    """
    mol = decode_mol(pmol)
    if mol is None:
        raise ToolError("Failed to decode the pickled RDKit Mol object.")

    pdb_string = Chem.MolToPDBBlock(mol)

    if filename is None:
        filename = f"{Chem.MolToSmiles(mol)}.pdb"
    if not filename.endswith('.pdb'):
        filename += '.pdb'
    with tempfile.NamedTemporaryFile(suffix=".pdb") as temp_pdb_file:
        temp_pdb_file.write(pdb_string.encode('utf-8'))
        temp_pdb_file.flush()
        return encode_file_contents(temp_pdb_file.name, filename=filename)


@rdkit_tool()
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


@rdkit_tool()
def sdf_contents_to_mol(sdf_contents: str) -> PickledMol:
    """
    Converts the contents of an SDF file to a pickled RDKit Mol object.

    Args:
        sdf_contents: The contents of the SDF file.
    Returns:
        A base64 encoded pickled RDKit Mol object.
    """
    sdf_io = BytesIO(sdf_contents.encode('utf-8'))
    supplier = Chem.ForwardSDMolSupplier(sdf_io)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None:
        raise ToolError(f"Failed to read molecule from SDF contents.")
    encoded_mol = encode_mol(mol)
    return encoded_mol
