import logging

from rdkit import Chem
from rdkit.Chem import rdMMPA
from mcp.server.fastmcp.exceptions import ToolError
from typing import List
from ..decorators import rdkit_tool
from ..types import Smiles, MolFragments

logger = logging.getLogger(__name__)


@rdkit_tool()
def FragmentMol(
    smiles: Smiles,
    maxCuts: int = 3,
    maxCutBonds: int = 20,
    pattern: str = "[#6+0;!$(*=,#[!#6])]!@!=!#[*]",
) -> List[MolFragments]:
    """Does the fragmentation necessary for an MMPA analysis.

    Parameters:
    - smiles (str): The SMILES string of the molecule to be fragmented.
    - maxCuts (int): Maximum number of cuts to make in the molecule.
    - maxCutBonds (int): Maximum number of bonds to cut.
    - pattern (str): An rSMARTS string defining the bond-breaking SMARTS pattern to use.

    Returns:
    - list[MolFragments]: A list of fragment pairs. Each pair is a 2-tuple
      `(core_smiles, side_chain_smiles)`, where:
        - `core_smiles` is either a SMILES string of the common scaffold or `None`
          if no core could be defined,
        - `side_chain_smiles` is always a valid SMILES string containing one or
          more wildcard atoms (e.g., `[*:1]`) indicating attachment points.

    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ToolError("Invalid SMILES string")

    pattern_mol = Chem.MolFromSmarts(pattern)
    if not pattern_mol:
        raise ToolError("Invalid pattern SMARTS")

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
    for core, side in frags:
        inner_tuple = (
            Chem.MolToSmiles(core) if isinstance(core, Chem.Mol) else None,
            Chem.MolToSmiles(side),
        )
        output.append(inner_tuple)
    logger.debug(f"Fragmentation output: {output}")
    return output
