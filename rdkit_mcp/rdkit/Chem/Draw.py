import os
from mcp.server.fastmcp.exceptions import ToolError

from rdkit import Chem
from rdkit.Chem import Draw


from ...tools.utils import rdkit_tool
from ...tools.types import Smiles


import logging
from rdkit.Chem.Draw import *

logger = logging.getLogger(__name__)


@rdkit_tool(description=Draw.MolToFile.__doc__)
def MolToFile(smiles: Smiles, filepath: str, width: int = 300, height: int = 300):
    mol: Chem.Mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ToolError(f"Invalid or unparsable SMILES string: {smiles}")

    if not filepath.endswith('.png'):
        filepath += '.png'

    output_path = os.path.join(filepath)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    Draw.MolToFile(mol, output_path, size=(width, height))
    return output_path
