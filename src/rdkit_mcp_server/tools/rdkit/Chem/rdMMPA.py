import logging

from rdkit import Chem
from rdkit.Chem import rdMMPA
from mcp.server.fastmcp.exceptions import ToolError
from typing import List
from ...utils import rdkit_tool
logger = logging.getLogger(__name__)


@rdkit_tool()
def FragmentMol(
    smiles: str,
    maxCuts: int = 3,
    maxCutBonds: int = 20,
    pattern: str = '',
) -> List[List[str]]:
    """Does the fragmentation necessary for an MMPA analysis.
    
    Parameters:
    - smiles (str): The SMILES string of the molecule to be fragmented.
    - maxCuts (int): Maximum number of cuts to make in the molecule.
    - maxCutBonds (int): Maximum number of bonds to cut.
    - pattern (str): The SMARTS pattern to use for fragmentation. Default is '[#6+0;!$(*=, #[!#6])]!@!=!#[*]'.
    
    Returns:
    - tuple: A tuple containing the fragments as RDKit Mol objects.
    
    """
    pattern = pattern or '[#6+0,#7+0;R;!$(*=,#[!#6])]!@!=!#[*]'
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ToolError("Invalid SMILES string")

    frags = rdMMPA.FragmentMol(
        mol,
        maxCuts=maxCuts,
        maxCutBonds=maxCutBonds,
        pattern=pattern,
        resultsAsMols=True,
    )
    if not frags:
        raise ToolError("Fragmentation failed")
    logger.info(f"Fragmentation produced {len(frags)} fragments")
    logger.info(frags)

    output = []
    for frag_tup in frags:
        inner_list = []
        for item in frag_tup:
            if isinstance(item, Chem.Mol):
                inner_list.append(Chem.MolToSmiles(item))
            else:
                inner_list.append(item)
        output.append(inner_list)
    return output
