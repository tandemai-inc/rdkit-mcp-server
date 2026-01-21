from typing import Annotated, Literal
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
    CalcNumSaturatedCarbocycles as _CalcNumSaturatedCarbocycles,
    CalcNumSaturatedRings as _CalcNumSaturatedRings,
    CalcNumSpiroAtoms as _CalcNumSpiroAtoms,
    CalcNumUnspecifiedAtomStereoCenters as _CalcNumUnspecifiedAtomStereoCenters,
    CalcOxidationNumbers as _CalcOxidationNumbers,
)

from ..decorators import rdkit_tool
from ..types import Smiles


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


@rdkit_tool(description=_CalcOxidationNumbers.__doc__)
def CalcOxidationNumbers(smiles: Smiles) -> float:
    """Calculate the oxidation numbers for a molecule given its SMILES representation."""
    mol = Chem.MolFromSmiles(smiles)
    return _CalcOxidationNumbers(mol)


# Names gathered from rdkit.Chem.rdMolDescriptors.Properties()
DescriptorNames = Annotated[
    Annotated[Literal["exactmw"], Field(description="Exact molecular weight (monoisotopic mass)")]
    | Annotated[Literal["amw"], Field(description="Average molecular weight")]
    | Annotated[Literal["lipinskiHBA"], Field(description="Hydrogen bond acceptors (Lipinski definition: N + O count)")]
    | Annotated[Literal["lipinskiHBD"], Field(description="Hydrogen bond donors (Lipinski definition: NH + OH count)")]
    | Annotated[Literal["NumRotatableBonds"], Field(description="Number of rotatable bonds (excludes terminal bonds)")]
    | Annotated[Literal["NumHBD"], Field(description="Number of hydrogen bond donors")]
    | Annotated[Literal["NumHBA"], Field(description="Number of hydrogen bond acceptors")]
    | Annotated[Literal["NumHeavyAtoms"], Field(description="Number of non-hydrogen atoms")]
    | Annotated[Literal["NumAtoms"], Field(description="Total number of atoms including hydrogens")]
    | Annotated[Literal["NumHeteroatoms"], Field(description="Number of non-carbon, non-hydrogen atoms")]
    | Annotated[Literal["NumAmideBonds"], Field(description="Number of amide bonds (C(=O)-N)")]
    | Annotated[Literal["FractionCSP3"], Field(description="Fraction of sp3 hybridized carbons")]
    | Annotated[Literal["NumRings"], Field(description="Total number of rings (SSSR count)")]
    | Annotated[Literal["NumAromaticRings"], Field(description="Number of aromatic rings")]
    | Annotated[Literal["NumAliphaticRings"], Field(description="Number of aliphatic (non-aromatic) rings")]
    | Annotated[Literal["NumSaturatedRings"], Field(description="Number of saturated rings (no double bonds)")]
    | Annotated[Literal["NumHeterocycles"], Field(description="Number of rings containing heteroatoms")]
    | Annotated[Literal["NumAromaticHeterocycles"], Field(description="Number of aromatic rings with heteroatoms")]
    | Annotated[Literal["NumSaturatedHeterocycles"], Field(description="Number of saturated rings with heteroatoms")]
    | Annotated[Literal["NumAliphaticHeterocycles"], Field(description="Number of aliphatic rings with heteroatoms")]
    | Annotated[Literal["NumSpiroAtoms"], Field(description="Number of spiro atoms (shared between two rings)")]
    | Annotated[Literal["NumBridgeheadAtoms"], Field(description="Number of bridgehead atoms in fused ring systems")]
    | Annotated[Literal["NumAtomStereoCenters"], Field(description="Number of atomic stereocenters (chiral centers)")]
    | Annotated[Literal["NumUnspecifiedAtomStereoCenters"], Field(description="Number of stereocenters with undefined stereochemistry")]
    | Annotated[Literal["labuteASA"], Field(description="Labute's approximate surface area")]
    | Annotated[Literal["tpsa"], Field(description="Topological polar surface area")]
    | Annotated[Literal["CrippenClogP"], Field(description="Wildman-Crippen LogP (lipophilicity estimate)")]
    | Annotated[Literal["CrippenMR"], Field(description="Wildman-Crippen molar refractivity")]
    | Annotated[Literal["chi0v"], Field(description="Chi0v molecular connectivity index (valence, order 0)")]
    | Annotated[Literal["chi1v"], Field(description="Chi1v molecular connectivity index (valence, order 1)")]
    | Annotated[Literal["chi2v"], Field(description="Chi2v molecular connectivity index (valence, order 2)")]
    | Annotated[Literal["chi3v"], Field(description="Chi3v molecular connectivity index (valence, order 3)")]
    | Annotated[Literal["chi4v"], Field(description="Chi4v molecular connectivity index (valence, order 4)")]
    | Annotated[Literal["kappa1"], Field(description="Kappa1 shape index (encodes atom count & branching)")]
    | Annotated[Literal["kappa2"], Field(description="Kappa2 shape index (encodes path/branching complexity)")]
    | Annotated[Literal["kappa3"], Field(description="Kappa3 shape index (encodes higher-order shape features)")],
    Field(description="Supported RDKit molecular property descriptor names."),
]


@rdkit_tool()
def compute_descriptors(
    smiles_list: Annotated[list[Smiles], Field(description="A list of SMILES strings representing ligands.")],
    descriptor_names: Annotated[list[DescriptorNames], Field(description="A list of RDKit molecular property names to compute")]
) -> list[list[float]]:
    """
    Compute RDKit molecular property descriptors for a collection of ligands.

    Each input SMILES string is parsed into an RDKit Mol object and all requested
    descriptors are computed in a single C++ call per molecule using
    ``rdMolDescriptors.Properties`` for efficiency.

    Parameters
    ----------
    smiles_list : list of str
        A list of SMILES strings representing ligands.

        - Each element must be a valid SMILES string accepted by
          ``Chem.MolFromSmiles``.
        - Invalid SMILES strings are silently skipped and do not appear in
          the output.

    descriptor_names : list of str
        A list of RDKit molecular property names to compute.

        - Each name must be supported by ``rdMolDescriptors.Properties``.
        - Property names are case-sensitive.
        - All descriptors are returned as scalar numeric values (floats).
        - Available property names include: exactmw, amw, NumHeavyAtoms, NumHBD,
          NumHBA, NumRotatableBonds, NumRings, NumAromaticRings, NumAliphaticRings,
          tpsa, CrippenClogP, CrippenMR, and more.

    Returns
    -------
    rows : list of list of float
        A list of descriptor rows, one per successfully parsed ligand.

        - Each inner list has length ``len(descriptor_names)``.
        - Descriptor values are returned in the same order as ``descriptor_names``.
        - Output order corresponds to the input order of valid SMILES only
          (invalid entries are omitted).

    Notes
    -----
    - Each ligand is parsed and sanitised exactly once.
    - All descriptors for a ligand are computed in a single call to minimise
      Python overhead.
    - This function does not perform multiprocessing, caching, or NaN filling.
    - No NumPy arrays are created; output is returned as plain Python lists.

    Raises
    ------
    ValueError
        If ``descriptor_names`` contains a name not supported by
        ``Chem.rdMolDescriptors.Properties``.
    """

    props = Chem.rdMolDescriptors.Properties(descriptor_names)

    rows = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        # Compute all descriptors in one call, convert to plain Python list
        rows.append(list(props.ComputeProperties(mol)))

    return rows
