from rdkit import Chem
from rdkit.Chem import Descriptors

from mcp.server.fastmcp.exceptions import ToolError

def exact_mol_wt(smiles: str) -> float:
    """
    Calculates the exact molecular weight of a molecule using its SMILES representation.

    Parameters:
    - smiles (str): The SMILES string of the molecule.

    Returns:
    - dict: A dictionary containing the exact molecular weight or an error message.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ToolError("Invalid SMILES string")
        weight = Descriptors.ExactMolWt(mol)
        return weight
    except Exception as e:
        raise ToolError(f"Error calculating ExactMolWt: {str(e)}")


def fp_density_morgan1(smiles: str) -> str:
    """
    Calculates the fingerprint density using Morgan fingerprints of radius 1.

    Parameters:
    - smiles (str): The SMILES string of the molecule.

    Returns:
    - float: A float containing the fingerprint density or an error message.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ToolError("Invalid SMILES string")
        density = Descriptors.FpDensityMorgan1(mol)
        return density
    except Exception as e:
        raise ToolError(f"Error calculating FpDensityMorgan1: {str(e)}")


def fp_density_morgan2(smiles: str) -> float:
    """
    Calculates the fingerprint density using Morgan fingerprints of radius 2.

    Parameters:
    - smiles (str): The SMILES string of the molecule.

    Returns:
    - float: A float containing the fingerprint density or an error message.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ToolError("Invalid SMILES string")
        density = Descriptors.FpDensityMorgan3(mol)
        return density
    except Exception as e:
         raise ToolError(f"Error calculating FpDensityMorgan3: {str(e)}")


def fp_density_morgan3(smiles: str) -> float:
    """
    Calculates the fingerprint density using Morgan fingerprints of radius 3.

    Parameters:
    - smiles (str): The SMILES string of the molecule.

    Returns:
    - float: A float containing the fingerprint density or an error message.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ToolError("Invalid SMILES string")
        density = Descriptors.FpDensityMorgan3(mol)
        return {"FpDensityMorgan3": density}
    except Exception as e:
         raise ToolError(f"Error calculating FpDensityMorgan3: {str(e)}")


def heavy_atom_mol_wt(smiles: str) -> float:
    """
    Calculates the molecular weight of a molecule excluding hydrogen atoms.

    Parameters:
    - smiles (str): The SMILES string of the molecule.

    Returns:
    - dict: A dictionary containing the heavy atom molecular weight or an error message.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ToolError("Invalid SMILES string")
        weight = Descriptors.HeavyAtomMolWt(mol)
        return {"HeavyAtomMolWt": weight}
    except Exception as e:
         raise ToolError(f"Error calculating HeavyAtomMolWt: {str(e)}")

def max_abs_partial_charge(smiles: str) -> dict:
    """
    Calculates the maximum absolute partial charge of a molecule.

    Parameters:
    - smiles (str): The SMILES string of the molecule.

    Returns:
    - dict: A dictionary containing the maximum absolute partial charge or an error message.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ToolError("Invalid SMILES string")
        charge = Descriptors.MaxAbsPartialCharge(mol)
        return {"MaxAbsPartialCharge": charge}
    except Exception as e:
         raise ToolError(f"Error calculating MaxAbsPartialCharge: {str(e)}")


def max_partial_charge(smiles: str) -> float:
    """
    Calculates the maximum partial charge of a molecule.

    Parameters:
    - smiles (str): The SMILES string of the molecule.

    Returns:
    - dict: A dictionary containing the maximum partial charge or an error message.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ToolError("Invalid SMILES string")
        charge = Descriptors.MaxPartialCharge(mol)
        return {"MaxPartialCharge": charge}
    except Exception as e:
         raise ToolError(f"Error calculating MaxPartialCharge: {str(e)}")


def min_abs_partial_charge(smiles: str) -> float:
    """
    Calculates the minimum absolute partial charge of a molecule.

    Parameters:
    - smiles (str): The SMILES string of the molecule.

    Returns:
    - dict: A dictionary containing the minimum absolute partial charge or an error message.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ToolError("Invalid SMILES string")
        charge = Descriptors.MinAbsPartialCharge(mol)
        return charge
    except Exception as e:
         raise ToolError(f"Error calculating MinAbsPartialCharge: {str(e)}")


def min_partial_charge(smiles: str) -> float:
    """
    Calculates the minimum partial charge of a molecule.

    Parameters:
    - smiles (str): The SMILES string of the molecule.

    Returns:
    - dict: A dictionary containing the minimum partial charge or an error message.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ToolError("Invalid SMILES string")
        charge = Descriptors.MinPartialCharge(mol)
        return charge
    except Exception as e:
        raise ToolError(f"Error calculating MinPartialCharge: {str(e)}")


def mol_wt(smiles: str) -> float:
    """
    Calculates the molecular weight of a molecule including hydrogen atoms.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ToolError("Invalid SMILES string")
        mol_wt =  Descriptors.MolWt(mol)
        return mol_wt
    except Exception as e:
        raise ToolError(str(e))


def num_radical_electrons(smiles: str) -> int:
    """
    Calculates the number of radical electrons in a molecule.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ToolError("Invalid SMILES string")
        return Descriptors.NumRadicalElectrons(mol)
    except Exception as e:
         raise ToolError(str(e))


def num_valence_electrons(smiles: str) -> int:
    """
    Calculates the number of valence electrons in a molecule.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
             raise ToolError("Invalid SMILES string")
        return Descriptors.NumValenceElectrons(mol)
    except Exception as e:
         raise ToolError(str(e))
