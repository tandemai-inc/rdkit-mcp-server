import logging
from typing import Tuple


from rdkit.Chem import Mol

from mcp.server.fastmcp.exceptions import ToolError
from ...decorators import rdkit_tool
from ...types import PickledMol, Smiles
from ...utils import encode_mol, decode_mol

logger = logging.getLogger(__name__)


@rdkit_tool(description=Mol.GetSubstructMatch.__doc__)
def GetSubstructMatch(
        p_mol: PickledMol,
        p_query: PickledMol) -> Tuple[int]:
    try:
        mol: Mol = decode_mol(p_mol)
        query: Mol = decode_mol(p_query)
        if not mol:
            raise ToolError("Invalid pickled RDKit Mol object")
        if not query:
            raise ToolError("Invalid pickled RDKit query Mol object")
        matching_atom_ids = mol.GetSubstructMatch(query)
        return tuple(matching_atom_ids)
    except Exception as e:
        raise ToolError(f"Error calculating GetSubstructMatch: {str(e)}")


@rdkit_tool(description=Mol.HasSubstructMatch.__doc__)
def HasSubstructMatch(
    p_mol: PickledMol,
    p_query: PickledMol,
    recursion_possible: bool = True,
    use_chirality: bool = False,
    use_query_query_matches: bool = False,
) -> bool:
    mol: Mol = decode_mol(p_mol)
    query: Mol = decode_mol(p_query)
    if not mol:
        raise ToolError("Invalid pickled RDKit Mol object")
    if not query:
        raise ToolError("Invalid pickled RDKit query Mol object")

    has_match: bool = mol.HasSubstructMatch(
        query,
        recursionPossible=recursion_possible,
        useChirality=use_chirality,
        useQueryQueryMatches=use_query_query_matches
    )
    return has_match


@rdkit_tool(description=Mol.SetProp.__doc__)
def SetProp(
    p_mol: PickledMol,
    key: str,
    value: str,
    computed: bool = False,
) -> PickledMol:
    mol: Mol = decode_mol(p_mol)
    mol.SetProp(key, value, computed)
    return encode_mol(mol)


@rdkit_tool(description=Mol.SetProp.__doc__)
def SetIntProp(
    p_mol: PickledMol,
    key: str,
    value: int,
    computed: bool = False,
) -> PickledMol:
    mol: Mol = decode_mol(p_mol)
    mol.SetIntProp(key, value, computed)
    return encode_mol(mol)


@rdkit_tool(description=Mol.SetProp.__doc__)
def SetBoolProp(
    p_mol: PickledMol,
    key: str,
    value: bool,
    computed: bool = False,
) -> PickledMol:
    mol: Mol = decode_mol(p_mol)
    mol.SetBoolProp(key, value, computed)
    return encode_mol(mol)


@rdkit_tool(description=Mol.SetProp.__doc__)
def SetDoubleProp(
    p_mol: PickledMol,
    key: str,
    value: float,
    computed: bool = False,
) -> PickledMol:
    mol: Mol = decode_mol(p_mol)
    mol.SetDoubleProp(key, value, computed)
    return encode_mol(mol)


@rdkit_tool(description=Mol.SetProp.__doc__)
def SetUnsignedProp(
    p_mol: PickledMol,
    key: str,
    value: int,
    computed: bool = False,
) -> PickledMol:
    mol: Mol = decode_mol(p_mol)
    mol.SetUnsignedProp(key, value, computed)
    return encode_mol(mol)


@rdkit_tool(description=Mol.SetProp.__doc__)
def UpdatePropertyCache(
    p_mol: PickledMol,
    key: str,
    value: int,
    computed: bool = False,
) -> PickledMol:
    mol: Mol = decode_mol(p_mol)
    mol.SetUnsignedProp(key, value, computed)
    return encode_mol(mol)
