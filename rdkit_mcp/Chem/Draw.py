from io import BytesIO
import logging
import os

from datetime import datetime
from mcp.server.fastmcp.exceptions import ToolError
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
from typing import List, Annotated

from ..decorators import rdkit_tool
from ..types import PickledMol, Smiles, EncodedFile
from ..utils import encode_file_contents
from rdkit_mcp.settings import ToolSettings
from rdkit_mcp.utils import decode_mol


logger = logging.getLogger(__name__)


@rdkit_tool(description=Draw.MolToFile.__doc__)
def MolToFile(pmol: PickledMol, filename: str, width: int = 300, height: int = 300) -> EncodedFile:
    mol: Chem.Mol = decode_mol(pmol)
    if mol is None:
        raise ToolError(f"Invalid or unparsable SMILES string: {smiles}")

    if not filename.endswith('.png'):
        filename += '.png'

    settings = ToolSettings()
    output_path = os.path.join(settings.FILE_DIR, filename)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    buffer = BytesIO()
    Draw.MolToFile(mol, buffer, size=(width, height))
    return encode_file_contents(buffer, filename=filename)


@rdkit_tool(description=Draw.MolsMatrixToGridImage.__doc__)
def MolsMatrixToGridImage(
    molsMatrix: List[List[Smiles]],
    subImgSize: Annotated[list[int], "2 dimensional image size for each sub-image in matrix"] = [200, 200],
    legendsMatrix: List[List[str]] = None,
    highlightAtomListsMatrix: List[List[int]] = None,
    highlightBondListsMatrix: List[List[int]] = None,
    useSVG: bool = False,
    returnPNG: bool = False,
    filename: Annotated[str, "output filename"] = None,
) -> EncodedFile:
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

    buffer = BytesIO()
    img.save(buffer, format="PNG")
    return encode_file_contents(buffer, filename=filename)


@rdkit_tool(description=Draw.MolToImage.__doc__)
def MolToImage(
    pmol: PickledMol,
    size: list[int, int] = [300, 300],
    kekulize: bool = True,
    wedgeBonds: bool = True,
    fitImage: bool = False,
    # options=None,
    filename: Annotated[str, "output filename"] = None,
    highlightAtoms: Annotated[list[int], "List of atom ids to highlight in image"] = None,
    highlightBonds: Annotated[list[int], "List of bond ids to highlight in image"] = None,
    highlightColor: Annotated[list[float], "Highlight RGB color"] = [1, 0, 0],
    **kwargs
) -> EncodedFile:
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
    mol: Chem.Mol = decode_mol(pmol)
    if mol is None:
        raise ToolError(f"Invalid or unparsable SMILES string: {smiles}")

    img = Draw.MolToImage(
        mol,
        size=size,
        kekulize=kekulize,
        wedgeBonds=wedgeBonds,
        fitImage=fitImage,
        # options=options,
        highlightAtoms=highlightAtoms,
        highlightBonds=highlightBonds,
        highlightColor=highlightColor,
        **kwargs
    )
    buffer = BytesIO()
    img.save(buffer, format="PNG")
    return encode_file_contents(buffer, filename=filename)
