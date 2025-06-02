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

from ...utils import rdkit_tool

smiles_type = Annotated[str, Field(description="SMILES representation of a molecule")]


@rdkit_tool(enabled=False)
def calc_exact_mol_wt(*args, **kwargs):
    return CalcExactMolWt(*args, **kwargs)


@rdkit_tool(enabled=True, description=CalcMolFormula.__doc__)
def calc_mol_formula(smiles: smiles_type) -> str:
    mol = Chem.MolFromSmiles(smiles)
    return CalcMolFormula(mol)


@rdkit_tool(enabled=True, description=CalcNumRings.__doc__)
def calc_num_rings(
    smiles: smiles_type,
) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return CalcNumRings(mol)


@rdkit_tool(enabled=True, description=CalcNumAromaticRings.__doc__)
def calc_num_aromatic_rings(
    smiles: smiles_type,
) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return CalcNumAromaticRings(mol)


@rdkit_tool(enabled=False)
def calc_num_aliphatic_rings(*args, **kwargs):
    return CalcNumAliphaticRings(*args, **kwargs)


@rdkit_tool(enabled=True, description=CalcNumRotatableBonds.__doc__)
def calc_num_rotatable_bonds(
        smiles: smiles_type,
        strict: Annotated[bool, Field(description="Handles linkages between ring systems.")] = True) -> int:
    """Calculate the number of rotatable bonds in a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return CalcNumRotatableBonds(mol, strict)


@rdkit_tool(enabled=False)
def calc_tpsa(*args, **kwargs):
    return CalcTPSA(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_chi0v(*args, **kwargs):
    return CalcChi0v(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_chi1v(*args, **kwargs):
    return CalcChi1v(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_chi2v(*args, **kwargs):
    return CalcChi2v(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_chi3v(*args, **kwargs):
    return CalcChi3v(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_chi4v(*args, **kwargs):
    return CalcChi4v(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_kappa1(*args, **kwargs):
    return CalcKappa1(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_kappa2(*args, **kwargs):
    return CalcKappa2(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_kappa3(*args, **kwargs):
    return CalcKappa3(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_labute_asa(*args, **kwargs):
    return CalcLabuteASA(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_crippen_descriptors(*args, **kwargs):
    return CalcCrippenDescriptors(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_num_hba(*args, **kwargs):
    return CalcNumHBA(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_num_hbd(*args, **kwargs):
    return CalcNumHBD(*args, **kwargs)


@rdkit_tool(enabled=True, description=CalcNumLipinskiHBA.__doc__)
def calc_num_lipinski_hba(
    smiles: smiles_type
) -> int:
    """Calculate the number of hydrogen bond acceptors according to Lipinski's rule of five."""
    mol = Chem.MolFromSmiles(smiles)
    return CalcNumLipinskiHBA(mol)


@rdkit_tool(enabled=True, description=CalcNumLipinskiHBD.__doc__)
def calc_num_lipinski_hbd(smiles: smiles_type) -> int:
    """Calculate the number of hydrogen bond donors according to Lipinski's rule of five."""
    mol = Chem.MolFromSmiles(smiles)
    return CalcNumLipinskiHBD(mol)


@rdkit_tool(enabled=False)
def calc_num_amide_bonds(*args, **kwargs):
    return CalcNumAmideBonds(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_fraction_csp3(*args, **kwargs):
    return CalcFractionCSP3(*args, **kwargs)


@rdkit_tool(enabled=False)
def calc_pbf(*args, **kwargs):
    return CalcPBF(*args, **kwargs)


@rdkit_tool(enabled=False)
def get_usr_score(*args, **kwargs):
    return GetUSRScore(*args, **kwargs)


@rdkit_tool(enabled=False)
def get_usr(*args, **kwargs):
    return GetUSR(*args, **kwargs)
