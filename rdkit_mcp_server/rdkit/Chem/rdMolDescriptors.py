from typing import Annotated
from pydantic import Field
from rdkit import Chem

from rdkit.Chem.rdMolDescriptors import (
    CalcExactMolWt as _CalcExactMolWt,
    CalcMolFormula as _CalcMolFormula,
    CalcNumRings as _CalcNumRings,
    CalcNumAromaticRings as _CalcNumAromaticRings,
    CalcNumAliphaticRings as _CalcNumAliphaticRings,
    CalcNumRotatableBonds as _CalcNumRotatableBonds,
    CalcTPSA as _CalcTPSA,
    CalcChi0v as _CalcChi0v,
    CalcChi1v as _CalcChi1v,
    CalcChi2v as _CalcChi2v,
    CalcChi3v as _CalcChi3v,
    CalcChi4v as _CalcChi4v,
    CalcKappa1 as _CalcKappa1,
    CalcKappa2 as _CalcKappa2,
    CalcKappa3 as _CalcKappa3,
    CalcLabuteASA as _CalcLabuteASA,
    CalcCrippenDescriptors as _CalcCrippenDescriptors,
    CalcNumHBA as _CalcNumHBA,
    CalcNumHBD as _CalcNumHBD,
    CalcNumLipinskiHBA as _CalcNumLipinskiHBA,
    CalcNumLipinskiHBD as _CalcNumLipinskiHBD,
    CalcNumAmideBonds as _CalcNumAmideBonds,
    CalcFractionCSP3 as _CalcFractionCSP3,
    CalcPBF as _CalcPBF,
    GetUSRScore as _GetUSRScore,
    GetUSR as _GetUSR,
    CalcNumHeterocycles as _CalcNumHeterocycles,
    CalcNumHeavyAtoms as _CalcNumHeavyAtoms,
    CalcNumHeteroatoms as _CalcNumHeteroatoms,
    CalcNumHeterocycles as _CalcNumHeterocycles,
    CalcNumSaturatedCarbocycles as _CalcNumSaturatedCarbocycles,
    CalcNumSaturatedRings as _CalcNumSaturatedRings,
    CalcNumSpiroAtoms as _CalcNumSpiroAtoms,
    CalcNumUnspecifiedAtomStereoCenters as _CalcNumUnspecifiedAtomStereoCenters,
    CalcOxidationNumbers as _CalcOxidationNumbers,
    CalcPBF as _CalcPBF,

)

from ...tools.utils import rdkit_tool
from ...tools.types import Smiles


@rdkit_tool(enabled=False)
def CalcExactMolWt(*args, **kwargs):
    return _CalcExactMolWt(*args, **kwargs)


@rdkit_tool(description=_CalcMolFormula.__doc__)
def CalcMolFormula(smiles: Smiles) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return _CalcMolFormula(mol)


@rdkit_tool(description=_CalcNumRings.__doc__)
def CalcNumRings(
    smiles: Smiles,
) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumRings(mol)


@rdkit_tool(description=_CalcNumAromaticRings.__doc__)
def CalcNumAromaticRings(
    smiles: Smiles,
) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumAromaticRings(mol)


@rdkit_tool()
def CalcNumAliphaticRings(smiles: Smiles) -> int:
    """Calculate the number of aliphatic rings in a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumAliphaticRings(mol)


@rdkit_tool(description=_CalcNumRotatableBonds.__doc__)
def CalcNumRotatableBonds(
        smiles: Smiles,
        strict: Annotated[bool, Field(description="Handles linkages between ring systems.")] = True) -> int:
    """Calculate the number of rotatable bonds in a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumRotatableBonds(mol, strict)


@rdkit_tool(description=_CalcTPSA.__doc__)
def CalcTPSA(smiles: Smiles) -> float:
    mol = Chem.MolFromSmiles(smiles)
    return _CalcTPSA(mol)


@rdkit_tool(description=_CalcChi0v.__doc__)
def CalcChi0v(smiles: Smiles) -> float:
    mol = Chem.MolFromSmiles(smiles)
    return _CalcChi0v(mol)


@rdkit_tool(description=_CalcChi1v.__doc__)
def CalcChi1v(smiles: Smiles) -> float:
    """Calculate the Chi1 value for a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcChi1v(mol)


@rdkit_tool(description=_CalcChi2v.__doc__)
def CalcChi2v(smiles: Smiles) -> float:
    mol = Chem.MolFromSmiles(smiles)
    return _CalcChi2v(mol)


@rdkit_tool(description=_CalcChi3v.__doc__)
def CalcChi3v(smiles: Smiles) -> float:
    mol = Chem.MolFromSmiles(smiles)
    return _CalcChi3v(mol)


@rdkit_tool(description=_CalcChi4v.__doc__)
def CalcChi4v(smiles: Smiles) -> float:
    """Calculate the Chi4 value for a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcChi4v(mol)


@rdkit_tool(description=_CalcKappa1.__doc__)
def CalcKappa1(smiles: Smiles) -> float:
    mol = Chem.MolFromSmiles(smiles)
    return _CalcKappa1(mol)


@rdkit_tool(description=_CalcKappa2.__doc__)
def CalcKappa2(smiles: Smiles) -> float:
    """Calculate the Kappa2 value for a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcKappa2(mol)


@rdkit_tool(description=_CalcKappa3.__doc__)
def CalcKappa3(smiles: Smiles) -> float:
    """Calculate the Kappa3 value for a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcKappa3(mol)


@rdkit_tool(description=_CalcLabuteASA.__doc__)
def CalcLabuteASA(smiles: Smiles, includeHs: bool = True, force: bool = False) -> float:
    """Calculate the Labute's solvent accessible surface area (ASA) for a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcLabuteASA(mol, includeHs=includeHs, force=force)


@rdkit_tool(description=_CalcCrippenDescriptors.__doc__)
def CalcCrippenDescriptors(smiles: Smiles, includeHs: bool = True, force: bool = False) -> tuple[float, float]:
    mol = Chem.MolFromSmiles(smiles)
    return _CalcCrippenDescriptors(mol, includeHs, force)


@rdkit_tool(description=_CalcNumHBA.__doc__)
def CalcNumHBA(smiles: Smiles) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumHBA(mol)


@rdkit_tool(description=_CalcNumHBD.__doc__)
def CalcNumHBD(smiles: Smiles) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumHBD(mol)


@rdkit_tool(description=_CalcNumLipinskiHBA.__doc__)
def CalcNumLipinskiHBA(
    smiles: Smiles
) -> int:
    """Calculate the number of hydrogen bond acceptors according to Lipinski's rule of five."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumLipinskiHBA(mol)


@rdkit_tool(description=_CalcNumLipinskiHBD.__doc__)
def CalcNumLipinskiHBD(smiles: Smiles) -> int:
    """Calculate the number of hydrogen bond donors according to Lipinski's rule of five."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumLipinskiHBD(mol)


@rdkit_tool(description=_CalcNumAmideBonds.__doc__)
def CalcNumAmideBonds(smiles: Smiles) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumAmideBonds(mol)


@rdkit_tool(description=_CalcNumHeterocycles.__doc__)
def CalcNumHeterocycles(smiles: Smiles) -> int:
    """Calculate the number of heterocycles in a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumHeterocycles(mol)


@rdkit_tool(description=_CalcNumHeavyAtoms.__doc__)
def CalcFractionCSP3(*args, **kwargs):
    return _CalcFractionCSP3(*args, **kwargs)


@rdkit_tool(description=_CalcPBF.__doc__)
def CalcPBF(*args, **kwargs):
    return _CalcPBF(*args, **kwargs)


@rdkit_tool(description=_GetUSRScore.__doc__)
def GetUSRScore(*args, **kwargs):
    return _GetUSRScore(*args, **kwargs)


@rdkit_tool(description=_GetUSR.__doc__)
def GetUSR(*args, **kwargs):
    return _GetUSR(*args, **kwargs)


@rdkit_tool(description=_CalcNumHeavyAtoms.__doc__)
def CalcNumHeavyAtoms(smiles: Smiles) -> int:
    """Calculate the number of heavy atoms in a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumHeavyAtoms(mol)


@rdkit_tool(description=_CalcNumHeteroatoms.__doc__)
def CalcNumHeteroatoms(smiles: Smiles) -> int:
    """Calculate the number of heteroatoms in a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumHeteroatoms(mol)


@rdkit_tool(description=_CalcNumHeterocycles.__doc__)
def CalcNumHeterocycles(smiles: Smiles) -> int:
    """Calculate the number of heterocycles in a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumHeterocycles(mol)


@rdkit_tool(description=_CalcNumSaturatedCarbocycles.__doc__)
def CalcNumSaturatedCarbocycles(smiles: Smiles) -> int:
    """Calculate the number of saturated carbocycles in a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumSaturatedCarbocycles(mol)


@rdkit_tool(description=_CalcNumSaturatedRings.__doc__)
def CalcNumSaturatedRings(smiles: Smiles) -> int:
    """Calculate the number of saturated rings in a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumSaturatedRings(mol)


@rdkit_tool(description=_CalcNumSpiroAtoms.__doc__)
def CalcNumSpiroAtoms(smiles: Smiles) -> int:
    """Calculate the number of spiro atoms in a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumSpiroAtoms(mol)


@rdkit_tool(description=_CalcNumUnspecifiedAtomStereoCenters.__doc__)
def CalcNumUnspecifiedAtomStereoCenters(smiles: Smiles) -> int:
    """Calculate the number of unspecified atom stereo centers in a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcNumUnspecifiedAtomStereoCenters(mol)


@rdkit_tool(description=_CalcPBF.__doc__)
def CalcPBF(smiles: Smiles) -> float:
    """Calculate the PBF (Polarizability, Basicity, and Flexibility) descriptor for a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcPBF(mol)
