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


@rdkit_tool()
def Compute2DCoords(
    p_mol: PickledMol,
    canonOrient: bool = True,
    clearConfs: bool = True,
    coordMap: dict = None,
    nFlipsPerSample: int = 0,
    nSample: int = 0,
    sampleSeed: int = 0,
    permuteDeg4Nodes: bool = False,
    bondLength: float = -1.0,
    forceRDKit: bool = False,
    useRingTemplates: bool = False
) -> PickledMol:
    """
    Compute 2D coordinates for a molecule.

    The resulting coordinates are stored on each atom of the molecule

    ARGUMENTS:

        mol - the molecule of interest canonOrient - orient the molecule in a canonical way clearConfs - if true, all existing conformations on the molecule

            will be cleared

        coordMap - a dictionary mapping atom Ids -> Point2D objects

            with starting coordinates for atoms that should have their positions locked.
        nFlipsPerSample - number of rotatable bonds that are

            flipped at random at a time.

        nSample - Number of random samplings of rotatable bonds. sampleSeed - seed for the random sampling process. permuteDeg4Nodes - allow permutation of bonds at a degree 4

            node during the sampling process

        bondLength - change the default bond length for depiction forceRDKit - use RDKit to generate coordinates even if

            preferCoordGen is set to true

        useRingTemplates - use templates to generate coordinates of complex

            ring systems

    RETURNS:
        ID of the conformation added to the molecule
    """
    mol: Mol = decode_mol(p_mol)
    mol.Compute2DCoords(
        canonOrient=canonOrient,
        clearConfs=clearConfs,
        coordMap=coordMap,
        nFlipsPerSample=nFlipsPerSample,
        nSample=nSample,
        sampleSeed=sampleSeed,
        permuteDeg4Nodes=permuteDeg4Nodes,
        bondLength=bondLength,
        forceRDKit=forceRDKit,
        useRingTemplates=useRingTemplates
    )
    return encode_mol(mol)
