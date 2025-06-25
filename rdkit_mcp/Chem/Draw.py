import os
from mcp.server.fastmcp.exceptions import ToolError
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
from typing import List

from ..decorators import rdkit_tool
from app.settings import get_app_settings
from ..types import Smiles

import logging
from rdkit.Chem.Draw import *

logger = logging.getLogger(__name__)


@rdkit_tool(description=Draw.MolToFile.__doc__)
def MolToFile(smiles: Smiles, filename: str, width: int = 300, height: int = 300) -> Path:
    mol: Chem.Mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ToolError(f"Invalid or unparsable SMILES string: {smiles}")

    if not filename.endswith('.png'):
        filename += '.png'

    settings = get_app_settings()
    output_path = os.path.join(settings.file_dir, filename)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    Draw.MolToFile(mol, output_path, size=(width, height))
    return Path(output_path)


@rdkit_tool(description=Draw.MolsMatrixToGridImage.__doc__)
def MolsMatrixToGridImage(
    molsMatrix: List[List[Smiles]],
    subImgSize: tuple = (200, 200),
    legendsMatrix: List[List[str]] = None,
    highlightAtomListsMatrix: List[List[int]] = None,
    highlightBondListsMatrix: List[List[int]] = None,
    useSVG: bool = False,
    returnPNG: bool = False,
    filename: str = None,
) -> Path:
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

    # Save the image to the specified file path
    settings = get_app_settings()
    file_path = os.path.join(settings.file_dir, filename)
    img.save(file_path)
    return Path(file_path)
