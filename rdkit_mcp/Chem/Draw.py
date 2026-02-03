import logging

from datetime import datetime
from pathlib import Path
from mcp.server.fastmcp.exceptions import ToolError
from rdkit import Chem
from rdkit.Chem import Draw
from typing import List, Annotated, Union

from ..decorators import rdkit_tool
from ..types import PickledMol, Smiles
from rdkit_mcp.utils import decode_mol


logger = logging.getLogger(__name__)


def _get_mol_from_inputs(
    pmol: PickledMol = None,
    pdb_path: Union[str, Path] = None,
    sdf_path: Union[str, Path] = None,
) -> Chem.Mol:
    """Load a molecule from either a PickledMol, PDB file, or SDF file."""
    inputs_provided = sum(x is not None for x in [pmol, pdb_path, sdf_path])

    if inputs_provided == 0:
        raise ToolError("Must provide one of: pmol, pdb_path, or sdf_path")
    if inputs_provided > 1:
        raise ToolError("Provide only one of: pmol, pdb_path, or sdf_path")

    if pmol is not None:
        mol = decode_mol(pmol)
        if mol is None:
            raise ToolError("Failed to decode pickled molecule")
        return mol

    if pdb_path is not None:
        pdb_path = Path(pdb_path)
        if not pdb_path.exists():
            raise ToolError(f"PDB file does not exist: {pdb_path}")
        mol = Chem.MolFromPDBFile(str(pdb_path))
        if mol is None:
            raise ToolError(f"Failed to read molecule from PDB file: {pdb_path}")
        return mol

    if sdf_path is not None:
        sdf_path = Path(sdf_path)
        if not sdf_path.exists():
            raise ToolError(f"SDF file does not exist: {sdf_path}")
        suppl = Chem.SDMolSupplier(str(sdf_path))
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            raise ToolError(f"Failed to read molecule from SDF file: {sdf_path}")
        return mol


@rdkit_tool(description=Draw.MolToFile.__doc__)
def MolToFile(
    file_dir: Union[str, Path],
    filename: str,
    pmol: Annotated[PickledMol, "Use only if you already have a PickledMol from a previous tool call"] = None,
    pdb_path: Annotated[Union[str, Path], "Path to PDB file - use this instead of calling pdb_to_mol first"] = None,
    sdf_path: Annotated[Union[str, Path], "Path to SDF file - use this instead of calling sdf_to_mol first"] = None,
    width: int = 300,
    height: int = 300,
) -> str:
    mol = _get_mol_from_inputs(pmol=pmol, pdb_path=pdb_path, sdf_path=sdf_path)

    if not filename.endswith('.png'):
        filename += '.png'

    filepath = Path(file_dir) / filename
    Draw.MolToFile(mol, filepath, size=(width, height))
    return str(filepath)


@rdkit_tool(description=Draw.MolsMatrixToGridImage.__doc__)
def MolsMatrixToGridImage(
    molsMatrix: List[List[Smiles]],
    subImgSize: Annotated[list[int], "2 dimensional image size for each sub-image in matrix"] = [200, 200],
    legendsMatrix: List[List[str]] = None,
    highlightAtomListsMatrix: List[List[int]] = None,
    highlightBondListsMatrix: List[List[int]] = None,
    useSVG: bool = False,
    returnPNG: bool = False,
    file_dir: Union[str, Path] = None,
    filename: Annotated[str, "output filename"] = None,
) -> str:
    if filename is None:
        raise ToolError("File path must be specified.")

    # Convert all SMILES in molsMatrix to Chem.Mol objects
    mol_matrix = []
    for row in molsMatrix:
        mol_row = []
        for smi in row:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                raise ToolError(f"Invalid or unparsable SMILES string: {smi}")
            mol_row.append(mol)
        mol_matrix.append(mol_row)
    molsMatrix = mol_matrix

    # Generate the grid image
    img = Draw.MolsMatrixToGridImage(
        molsMatrix,
        subImgSize=subImgSize,
        legendsMatrix=legendsMatrix,
        highlightAtomListsMatrix=highlightAtomListsMatrix,
        highlightBondListsMatrix=highlightBondListsMatrix,
        useSVG=useSVG,
        returnPNG=returnPNG)

    filepath = Path(file_dir) / filename
    img.save(filepath, format="PNG")
    return str(filepath)


@rdkit_tool(description=Draw.MolToImage.__doc__)
def MolToImage(
    file_dir: Union[str, Path],
    pmol: Annotated[PickledMol, "Base64 encoded pickled mol (use for small molecules)"] = None,
    pdb_path: Annotated[Union[str, Path], "Path to PDB file (use for large proteins)"] = None,
    sdf_path: Annotated[Union[str, Path], "Path to SDF file"] = None,
    size: list[int, int] = [300, 300],
    kekulize: bool = True,
    wedgeBonds: bool = True,
    fitImage: bool = False,
    filename: Annotated[str, "output filename"] = None,
    highlightAtoms: Annotated[list[int], "List of atom ids to highlight in image"] = None,
    highlightBonds: Annotated[list[int], "List of bond ids to highlight in image"] = None,
    highlightColor: Annotated[list[float], "Highlight RGB color"] = [1, 0, 0],
) -> str:
    if highlightAtoms is None:
        highlightAtoms = ()
    if highlightBonds is None:
        highlightBonds = ()
    if highlightColor is None:
        highlightColor = (1, 0, 0)
    if isinstance(highlightAtoms, list):
        highlightAtoms = tuple(highlightAtoms)
    if isinstance(highlightBonds, list):
        highlightBonds = tuple(highlightBonds)
    if isinstance(highlightColor, list):
        highlightColor = tuple(highlightColor)

    if filename is None:
        filename = f"mol_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"

    mol = _get_mol_from_inputs(pmol=pmol, pdb_path=pdb_path, sdf_path=sdf_path)

    img = Draw.MolToImage(
        mol,
        size=size,
        kekulize=kekulize,
        wedgeBonds=wedgeBonds,
        fitImage=fitImage,
        highlightAtoms=highlightAtoms,
        highlightBonds=highlightBonds,
        highlightColor=highlightColor,
    )
    filepath = Path(file_dir) / filename
    img.save(filepath, format="PNG")
    return str(filepath)
