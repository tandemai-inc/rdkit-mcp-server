from typing import Annotated
from pydantic import Field
from rdkit import Chem

from rdkit.Chem.rdMolDescriptors import (
    CalcExactMolWt,
    CalcMolFormula,
    CalcNumRings,
    CalcNumAromaticRings,
    CalcNumAliphaticRings,
    CalcNumRotatableBonds,
    CalcTPSA,
    CalcChi0v,
    CalcChi1v,
    CalcChi2v,
    CalcChi3v,
    CalcChi4v,
    CalcKappa1,
    CalcKappa2,
    CalcKappa3,
    CalcLabuteASA,
    CalcCrippenDescriptors,
    CalcNumHBA,
    CalcNumHBD,
    CalcNumLipinskiHBA,
    CalcNumLipinskiHBD,
    CalcNumAmideBonds,
    CalcFractionCSP3,
    CalcPBF,
    GetUSRScore,
    GetUSR,
)

from ...tools.utils import rdkit_tool

smiles_type = Annotated[str, Field(description="SMILES representation of a molecule")]


@rdkit_tool(enabled=False)
def CalcExactMolWt(*args, **kwargs):
    return CalcExactMolWt(*args, **kwargs)


@rdkit_tool(enabled=True, description=CalcMolFormula.__doc__)
def CalcMolFormula(smiles: smiles_type) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return CalcMolFormula(mol)


@rdkit_tool(enabled=True, description=CalcNumRings.__doc__)
def CalcNumRings(
    smiles: smiles_type,
) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return CalcNumRings(mol)


@rdkit_tool(enabled=True, description=CalcNumAromaticRings.__doc__)
def CalcNumAromaticRings(
    smiles: smiles_type,
) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return CalcNumAromaticRings(mol)


@rdkit_tool(enabled=False)
def calc_num_aliphatic_rings(*args, **kwargs):
    return CalcNumAliphaticRings(*args, **kwargs)


@rdkit_tool(enabled=True, description=CalcNumRotatableBonds.__doc__)
def CalcNumRotatableBonds(
        smiles: smiles_type,
        strict: Annotated[bool, Field(description="Handles linkages between ring systems.")] = True) -> int:
    """Calculate the number of rotatable bonds in a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return CalcNumRotatableBonds(mol, strict)


@rdkit_tool(enabled=False)
def CalcTPSA(*args, **kwargs):
    return CalcTPSA(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcChi0v(*args, **kwargs):
    return CalcChi0v(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcChi1v(*args, **kwargs):
    return CalcChi1v(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcChi2v(*args, **kwargs):
    return CalcChi2v(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcChi3v(*args, **kwargs):
    return CalcChi3v(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcChi4v(*args, **kwargs):
    return CalcChi4v(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcKappa1(*args, **kwargs):
    return CalcKappa1(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcKappa2(*args, **kwargs):
    return CalcKappa2(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcKappa3(*args, **kwargs):
    return CalcKappa3(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcLabuteASA(*args, **kwargs):
    return CalcLabuteASA(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcCrippenDescriptors(*args, **kwargs):
    return CalcCrippenDescriptors(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_num_hba(*args, **kwargs):
    return CalcNumHBA(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcNumHBD(*args, **kwargs):
    return CalcNumHBD(*args, **kwargs)


@rdkit_tool(enabled=True, description=CalcNumLipinskiHBA.__doc__)
def CalcNumLipinskiHBA(
    smiles: smiles_type
) -> int:
    """Calculate the number of hydrogen bond acceptors according to Lipinski's rule of five."""
    mol = Chem.MolFromSmiles(smiles)
    return CalcNumLipinskiHBA(mol)


@rdkit_tool(enabled=True, description=CalcNumLipinskiHBD.__doc__)
def CalcNumLipinskiHBD(smiles: smiles_type) -> int:
    """Calculate the number of hydrogen bond donors according to Lipinski's rule of five."""
    mol = Chem.MolFromSmiles(smiles)
    return CalcNumLipinskiHBD(mol)


@rdkit_tool(enabled=False)
def CalcNumAmideBonds(*args, **kwargs):
    return CalcNumAmideBonds(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcFractionCSP3(*args, **kwargs):
    return CalcFractionCSP3(*args, **kwargs)


@rdkit_tool(enabled=False)
def CalcPBF(*args, **kwargs):
    return CalcPBF(*args, **kwargs)


@rdkit_tool(enabled=False)
def GetUSRScore(*args, **kwargs):
    return GetUSRScore(*args, **kwargs)


@rdkit_tool(enabled=False)
def GetUSR(*args, **kwargs):
    return GetUSR(*args, **kwargs)
