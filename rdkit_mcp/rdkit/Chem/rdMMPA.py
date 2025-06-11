import logging

from rdkit import Chem
from rdkit.Chem import rdMMPA
from mcp.server.fastmcp.exceptions import ToolError
from typing import List
from ...tools.utils import rdkit_tool
from ...tools.types import Smiles, MolFragments
logger = logging.getLogger(__name__)


@rdkit_tool()
def FragmentMol(
    smiles: Smiles,
    maxCuts: int = 3,
    maxCutBonds: int = 20,
    pattern: str = '[#6+0;!$(*=,#[!#6])]!@!=!#[*]'
) -> List[MolFragments]:
    """Does the fragmentation necessary for an MMPA analysis.

    Parameters:
    - smiles (str): The SMILES string of the molecule to be fragmented.
    - maxCuts (int): Maximum number of cuts to make in the molecule.
    - maxCutBonds (int): Maximum number of bonds to cut.
    - pattern (str): An rSMARTS string defining the bond-breaking SMARTS pattern to use.

    Returns:
    - tuple: A tuple containing the fragments as RDKit Mol objects.

    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ToolError("Invalid SMILES string")
    try:
        frags = rdMMPA.FragmentMol(
            mol,
            maxCuts=maxCuts,
            maxCutBonds=maxCutBonds,
            pattern=pattern,
            resultsAsMols=True,
        )
    except Exception as e:
        logger.error(f"Fragmentation failed: {e}")
        raise ToolError(f"Fragmentation failed: {e}")
    if not frags:
        raise ToolError("Fragmentation failed")
    logger.debug(f"Fragmentation produced {len(frags)} fragments")
    output = []
    for frag_tup in frags:
        inner_list = []
        for item in frag_tup:
            if isinstance(item, Chem.Mol):
                inner_list.append(Chem.MolToSmiles(item))
            else:
                inner_list.append(item)
        output.append(inner_list)
    logger.debug(f"Fragmentation output: {output}")
    return output
