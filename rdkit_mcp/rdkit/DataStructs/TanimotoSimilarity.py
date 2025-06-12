import asyncio
from typing import Dict, Union, Optional
from mcp.server.fastmcp.exceptions import ToolError
from rdkit.Chem import AllChem
from rdkit import Chem, DataStructs

from ...tools.utils import rdkit_tool
from ...tools.types import Smiles
from rdkit.DataStructs import ExplicitBitVect


def compute_fingerprint(
    mol: Chem.Mol,
    radius: int = 2,
    nBits: int = 2048
) -> ExplicitBitVect:
    """
    Computes a molecular fingerprint for a given SMILES string.

    Args:
        smiles: The SMILES representation of the molecule.
        radius: The radius for Morgan fingerprints (default: 2). Ignored for 'rdkit'.
        nBits: The number of bits for the fingerprint (default: 2048).

    Returns:
        ExplicitBitVect: The computed fingerprint as an ExplicitBitVect object.
    """
    try:
        # Compute fingerprint
        fp: ExplicitBitVect = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
        return fp
    except Exception as e:
        raise ToolError(f"Error computing fingerprint: {e}")


@rdkit_tool()
async def TanimotoSimilarity(smiles1: Smiles, smiles2: Smiles, radius: int = 2, nBits: int = 2048) -> float:
    """
    Calculates the Tanimoto similarity between two molecules based on their fingerprints.

    Args:
        smiles1: The SMILES representation of the first molecule.
        smiles2: The SMILES representation of the second molecule.
        radius: Morgan fingerprint radius (default: 2).
        nBits: Fingerprint size in bits (default: 2048).

    Returns:
        A dictionary containing 'similarity_score' (a float between 0.0 and 1.0)
        and 'method' used, or an 'error' message string if calculation fails.
    """
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None:
        raise ToolError(f"Invalid or unparsable SMILES string for molecule 1: {smiles1}")
    if mol2 is None:
        raise ToolError(f"Invalid or unparsable SMILES string for molecule 2: {smiles2}")

    try:
        fp1: ExplicitBitVect = compute_fingerprint(mol1, radius, nBits)
        fp2: ExplicitBitVect = compute_fingerprint(mol2, radius, nBits)
        similarity = await asyncio.to_thread(DataStructs.TanimotoSimilarity, fp1, fp2)
        return round(similarity, 4)
    except ValueError as ve:
        raise ToolError(str(ve))
    except Exception as e:
        raise ToolError(f"Error calculating similarity: {e}")
