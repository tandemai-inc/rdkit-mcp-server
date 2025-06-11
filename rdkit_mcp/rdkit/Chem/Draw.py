import os
from mcp.server.fastmcp.exceptions import ToolError
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
from typing import List

from ...tools.utils import rdkit_tool
from ...tools.types import Smiles

import logging
from rdkit.Chem.Draw import *

logger = logging.getLogger(__name__)


@rdkit_tool(description=Draw.MolToFile.__doc__)
def MolToFile(smiles: Smiles, filepath: str, width: int = 300, height: int = 300) -> Path:
    mol: Chem.Mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ToolError(f"Invalid or unparsable SMILES string: {smiles}")

    if not filepath.endswith('.png'):
        filepath += '.png'

    output_path = os.path.join(filepath)
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
    file_path: str = None,
) -> Path:
    if file_path is None:
        raise ToolError("File path must be specified.")

    # Create the output directory if it doesn't exist
    if os.path.dirname(file_path):
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
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
    img.save(file_path)
    return Path(file_path)
