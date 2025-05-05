import asyncio
import logging
import os
import tempfile
import uuid
from typing import Dict, Union, Optional

# Import the MCP instance from the server module to use the @mcp.tool() decorator
from .server import mcp, logger

# Attempt to import RDKit and Pillow, logging errors if they are not found
try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem, Descriptors, Draw
    RDKIT_AVAILABLE = True
except ImportError:
    logger.error("RDKit library not found. Please install it (e.g., via conda). Tool functionality will be limited.")
    RDKIT_AVAILABLE = False

try:
    from PIL import Image
    PILLOW_AVAILABLE = True
except ImportError:
    logger.warning("Pillow library not found. Molecule drawing to PNG will not work.")
    PILLOW_AVAILABLE = False

# Helper function to handle RDKit molecule loading and errors
def _load_molecule(smiles: str) -> Optional[Chem.Mol]:
    """Loads a molecule from SMILES, returning None on failure."""
    if not RDKIT_AVAILABLE:
        logger.error("RDKit is not available.")
        return None
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

@mcp.tool()
async def parse_molecule(smiles: str) -> Dict[str, Union[str, int, float]]:
    """
    Parse a SMILES string into an RDKit molecule object and return basic properties.

    Args:
        smiles: The SMILES representation of the molecule.

    Returns:
        A dictionary containing 'atom_count', 'heavy_atom_count', 'formula',
        and 'molecular_weight' if parsing succeeds, or an 'error' message string if it fails.
    """
    logger.info(f"Tool 'parse_molecule' called with SMILES: {smiles[:30]}...")
    mol = await asyncio.to_thread(_load_molecule, smiles) # Run sync RDKit call in thread

    if mol is None:
        return {"error": f"Invalid or unparsable SMILES string: {smiles}"}

    try:
        # Calculate properties using RDKit (run in thread pool as they might block)
        atom_count = await asyncio.to_thread(mol.GetNumAtoms)
        heavy_atom_count = await asyncio.to_thread(mol.GetNumHeavyAtoms)
        # formula = await asyncio.to_thread(Descriptors.MolFormula, mol)
        mol_weight = await asyncio.to_thread(Descriptors.MolWt, mol)

        return {
            "atom_count": atom_count,
            "heavy_atom_count": heavy_atom_count,
            # "formula": formula,
            "molecular_weight": round(mol_weight, 4)
        }
    except Exception as e:
        logger.error(f"Error calculating properties for SMILES '{smiles}': {e}")
        return {"error": f"Error calculating properties: {e}"}

@mcp.tool()
async def draw_molecule(smiles: str, width: int = 300, height: int = 300) -> Dict[str, str]:
    """
    Generates a PNG image representation of a molecule from its SMILES string.

    Args:
        smiles: The SMILES representation of the molecule.
        width: The desired width of the image in pixels (default: 300).
        height: The desired height of the image in pixels (default: 300).

    Returns:
        A dictionary containing a 'file_uri' key with the file:// URI of the generated PNG image,
        or an 'error' key with an error message string if generation fails.
        Requires RDKit and Pillow libraries.
    """
    logger.info(f"Tool 'draw_molecule' called for SMILES: {smiles[:30]}...")
    if not RDKIT_AVAILABLE or not PILLOW_AVAILABLE:
        return {"error": "RDKit or Pillow library not available for drawing."}

    mol = await asyncio.to_thread(_load_molecule, smiles)
    if mol is None:
        return {"error": f"Invalid or unparsable SMILES string: {smiles}"}

    try:
        # Generate image using RDKit (sync call, run in thread)
        img = await asyncio.to_thread(Draw.MolToImage, mol, size=(width, height))

        # Save image to a temporary file
        temp_dir = tempfile.gettempdir()
        file_name = f"rdkit_mol_{uuid.uuid4()}.png"
        file_path = os.path.join(temp_dir, file_name)

        await asyncio.to_thread(img.save, file_path)

        # Generate file URI (ensure correct format for OS)
        # Note: Standard file URI format is file:///path/to/file
        # On Windows, it might be file:///C:/path/to/file
        if os.name == 'nt': # Windows
             file_uri = f"file:///{file_path.replace(os.sep, '/')}"
        else: # POSIX (macOS, Linux)
             file_uri = f"file://{file_path}"


        logger.info(f"Molecule image saved to: {file_path}")
        return {"file_uri": file_uri}

    except Exception as e:
        logger.error(f"Error drawing molecule for SMILES '{smiles}': {e}")
        return {"error": f"Error generating molecule image: {e}"}


@mcp.tool()
async def compute_fingerprint(smiles: str, method: str = "morgan", radius: int = 2, nBits: int = 2048) -> Dict[str, str]:
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
        return {"error": f"Invalid or unparsable SMILES string: {smiles}"}

    try:
        fp = None
        method_lower = method.lower()

        # Compute fingerprint (sync call, run in thread)
        if method_lower == "morgan":
            fp = await asyncio.to_thread(AllChem.GetMorganFingerprintAsBitVect, mol, radius, nBits=nBits)
        elif method_lower == "rdkit":
            fp = await asyncio.to_thread(Chem.RDKFingerprint, mol, fpSize=nBits)
        else:
            return {"error": f"Unsupported fingerprint method: {method}. Supported methods: 'morgan', 'rdkit'."}

        # Convert fingerprint to hex string
        fp_hex = await asyncio.to_thread(fp.ToBitString) # Get binary string first
        # Convert binary string to hex for better readability/storage if needed, though binary might be more standard
        # For simplicity, let's return the binary string representation directly.
        # fp_hex = fp.ToHex() # RDKit >= 2021.09 provides ToHex()
        # If using older RDKit or want hex from binary string:
        # fp_hex = hex(int(fp.ToBitString(), 2))[2:].upper() # Convert binary string to int, then hex

        return {"fingerprint_bits": fp.ToBitString(), "method": method_lower}

    except Exception as e:
        logger.error(f"Error computing fingerprint for SMILES '{smiles}': {e}")
        return {"error": f"Error computing fingerprint: {e}"}


@mcp.tool()
async def tanimoto_similarity(smiles1: str, smiles2: str, method: str = "morgan", radius: int = 2, nBits: int = 2048) -> Dict[str, Union[str, float]]:
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
        return {"error": f"Invalid or unparsable SMILES string for molecule 1: {smiles1}"}
    if mol2 is None:
        return {"error": f"Invalid or unparsable SMILES string for molecule 2: {smiles2}"}

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

    except ValueError as ve: # Catch specific error from get_fp
         logger.error(f"Fingerprint method error: {ve}")
         return {"error": str(ve)}
    except Exception as e:
        logger.error(f"Error calculating similarity between '{smiles1}' and '{smiles2}': {e}")
        return {"error": f"Error calculating similarity: {e}"}

logger.info("RDKit tool functions defined.")
