from io import BytesIO
import logging
from typing import Union
from mcp.server.fastmcp.exceptions import ToolError
from pathlib import Path
from rdkit import Chem

from .decorators import rdkit_tool
from .types import PickledMol, Smiles
from .utils import encode_mol, decode_mol

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
def mol_to_sdf(pmol: PickledMol, file_dir: Union[str, Path], filename: Union[str, None] = None) -> str:
    """
    Writes a pickled RDKit Mol object to an SDF file and returns its base64-encoded contents.

    Args:
        pmol: The pickled and base64-encoded RDKit Mol object.
        file_dir: Directory where the SDF file will be saved. If not provided, a temporary file is used.
        filename: Optional filename for the SDF file. If not provided, a name will be generated based on the SMILES string. Must end with .sdf
    Returns:
        Full filepath to the saved SDF File
    """
    mol = decode_mol(pmol)
    if mol is None:
        raise ToolError("Failed to decode the pickled RDKit Mol object.")

    if not file_dir:
        raise ToolError("file_dir must be provided and cannot be None.")

    file_dir = Path(file_dir)
    if filename is None:
        filename = f"{Chem.MolToSmiles(mol)}.sdf"
    if not filename.endswith('.sdf'):
        filename += '.sdf'
    file_path = file_dir / filename
    with Chem.SDWriter(str(file_path)) as writer:
        writer.write(mol)
    return str(file_path)


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
def mol_to_pdb(pmol: PickledMol, file_dir: Union[str, Path], filename: Union[str, None] = None) -> str:
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
    if not file_dir:
        raise ToolError("file_dir must be provided and cannot be None.")

    if filename is None:
        filename = f"{Chem.MolToSmiles(mol)}.pdb"
    if not filename.endswith('.pdb'):
        filename += '.pdb'

    filepath = Path(file_dir) / filename

    with Chem.PDBWriter(str(filepath)) as writer:
        writer.write(mol)
    return str(filepath)


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
