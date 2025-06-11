import asyncio
import os
from typing import Dict, Union, Optional
from mcp.server.fastmcp.exceptions import ToolError
from pathlib import Path
from rdkit.Chem import AllChem
from rdkit import Chem, DataStructs

from .utils import rdkit_tool, OUTPUT_DIR
from .types import Smiles
import logging

logger = logging.getLogger(__name__)

# Helper function to handle RDKit molecule loading and errors


def _load_molecule(smiles: Smiles) -> Optional[Chem.Mol]:
    """Loads a molecule from SMILES, returning None on failure."""
    if not smiles or not isinstance(smiles, str):
        logger.error("Invalid SMILES input: must be a non-empty string.")
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Failed to parse SMILES: {smiles}")
            return None
        return mol
    except Exception as e:
        logger.error(f"Error parsing SMILES '{smiles}': {e}")
        return None


@rdkit_tool()
async def compute_fingerprint(smiles: Smiles, method: str = "morgan", radius: int = 2, nBits: int = 2048) -> Dict[str, str]:
    """
    Computes a molecular fingerprint for a given SMILES string.

    Args:
        smiles: The SMILES representation of the molecule.
        method: The fingerprint method to use. Currently supports 'morgan' (default)
                or 'rdkit' (topological).
        radius: The radius for Morgan fingerprints (default: 2). Ignored for 'rdkit'.
        nBits: The number of bits for the fingerprint (default: 2048).

    Returns:
        A dictionary containing 'fingerprint_hex' (the fingerprint as a hex string)
        and 'method' used, or an 'error' message string if computation fails.
    """
    logger.info(f"Tool 'compute_fingerprint' called for SMILES: {smiles[:30]}..., Method: {method}")
    mol = await asyncio.to_thread(_load_molecule, smiles)
    if mol is None:
        raise ToolError(f"Invalid or unparsable SMILES string: {smiles}")

    try:
        fp = None
        method_lower = method.lower()

        # Compute fingerprint (sync call, run in thread)
        if method_lower == "morgan":
            fp = await asyncio.to_thread(AllChem.GetMorganFingerprintAsBitVect, mol, radius, nBits=nBits)
        elif method_lower == "rdkit":
            fp = await asyncio.to_thread(Chem.RDKFingerprint, mol, fpSize=nBits)
        else:
            raise ToolError(f"Unsupported fingerprint method: {method}. Supported methods: 'morgan', 'rdkit'.")

        # Convert fingerprint to hex string
        # fp_hex = await asyncio.to_thread(fp.ToBitString)  # Get binary string first
        # Convert binary string to hex for better readability/storage if needed, though binary might be more standard
        # For simplicity, let's return the binary string representation directly.
        # fp_hex = fp.ToHex() # RDKit >= 2021.09 provides ToHex()
        # If using older RDKit or want hex from binary string:
        # fp_hex = hex(int(fp.ToBitString(), 2))[2:].upper() # Convert binary string to int, then hex

        return {"fingerprint_bits": fp.ToBitString(), "method": method_lower}

    except Exception as e:
        logger.error(f"Error computing fingerprint for SMILES '{smiles}': {e}")
        raise ToolError(f"Error computing fingerprint: {e}")


@rdkit_tool()
async def tanimoto_similarity(smiles1: Smiles, smiles2: Smiles, method: str = "morgan", radius: int = 2, nBits: int = 2048) -> Dict[str, Union[str, float]]:
    """
    Calculates the Tanimoto similarity between two molecules based on their fingerprints.

    Args:
        smiles1: The SMILES representation of the first molecule.
        smiles2: The SMILES representation of the second molecule.
        method: The fingerprint method ('morgan' or 'rdkit') to use for comparison (default: 'morgan').
        radius: Morgan fingerprint radius (default: 2).
        nBits: Fingerprint size in bits (default: 2048).

    Returns:
        A dictionary containing 'similarity_score' (a float between 0.0 and 1.0)
        and 'method' used, or an 'error' message string if calculation fails.
    """
    logger.info(f"Tool 'tanimoto_similarity' called for SMILES1: {smiles1[:30]}..., SMILES2: {smiles2[:30]}..., Method: {method}")

    # Load molecules concurrently
    mol1_task = asyncio.to_thread(_load_molecule, smiles1)
    mol2_task = asyncio.to_thread(_load_molecule, smiles2)
    mol1, mol2 = await asyncio.gather(mol1_task, mol2_task)

    if mol1 is None:
        raise ToolError(f"Invalid or unparsable SMILES string for molecule 1: {smiles1}")
    if mol2 is None:
        raise ToolError(f"Invalid or unparsable SMILES string for molecule 2: {smiles2}")

    try:
        fp1 = None
        fp2 = None
        method_lower = method.lower()

        # Define fingerprinting function to run in thread
        def get_fp(mol, method, radius, nBits):
            if method == "morgan":
                return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
            elif method == "rdkit":
                return Chem.RDKFingerprint(mol, fpSize=nBits)
            else:
                raise ValueError(f"Unsupported fingerprint method: {method}")

        # Compute fingerprints concurrently
        fp1_task = asyncio.to_thread(get_fp, mol1, method_lower, radius, nBits)
        fp2_task = asyncio.to_thread(get_fp, mol2, method_lower, radius, nBits)
        fp1, fp2 = await asyncio.gather(fp1_task, fp2_task)

        # Calculate Tanimoto similarity (sync call, run in thread)
        similarity = await asyncio.to_thread(DataStructs.TanimotoSimilarity, fp1, fp2)

        return {"similarity_score": round(similarity, 4), "method": method_lower}

    except ValueError as ve:  # Catch specific error from get_fp
        logger.error(f"Fingerprint method error: {ve}")
        raise ToolError(str(ve))
    except Exception as e:
        logger.error(f"Error calculating similarity between '{smiles1}' and '{smiles2}': {e}")
        raise ToolError(f"Error calculating similarity: {e}")


@rdkit_tool()
def smiles_to_sdf(smiles: Smiles) -> Path:
    """
    Converts a SMILES string to an SDF file.

    Args:
        smiles: The SMILES representation of the molecule.

    Returns:
        An SDF string representation of the molecule.
    """
    logger.info(f"Converting SMILES to SDF: {smiles[:30]}...")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ToolError(f"Invalid or unparsable SMILES string: {smiles}")

    sdf_string = Chem.MolToMolBlock(mol)
    # Write to SDF file
    filename = f"{Chem.MolToSmiles(mol)}.sdf"
    output_path = Path(os.path.join(OUTPUT_DIR, filename))
    with open(output_path, "w") as f:
        f.write(sdf_string)
    return output_path


@rdkit_tool()
def sdf_to_smiles(sdf_path: Union[str, Path]) -> Smiles:
    """
    Converts an SDF file to a SMILES string.

    Args:
        sdf_path: The path to the SDF file.

    Returns:
        The SMILES representation of the molecule.
    """
    logger.info(f"Converting SDF to SMILES: {sdf_path}")
    if isinstance(sdf_path, str):
        sdf_path = Path(sdf_path)

    if not sdf_path.exists():
        raise ToolError(f"SDF file does not exist: {sdf_path}")

    suppl = Chem.SDMolSupplier(str(sdf_path))
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ToolError(f"Failed to read any valid molecule from SDF file: {sdf_path}")

    smiles = Chem.MolToSmiles(mol)
    return smiles


def get_base_tools():
    """
    Get all base tools defined in this module.

    Returns:
        A list of tool functions.
    """
    return [
        compute_fingerprint,
        tanimoto_similarity,
        smiles_to_sdf,
        sdf_to_smiles,
    ]
