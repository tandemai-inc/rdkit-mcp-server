import logging

from rdkit.Chem import Mol
from rdkit.Chem import rdDepictor as _rdDepictor
from ..decorators import rdkit_tool
from ..types import Smiles, PickledMol
from ..utils import encode_mol, decode_mol

logger = logging.getLogger(__name__)


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
    coordMap = coordMap or {}
    mol: Mol = decode_mol(p_mol)
    _rdDepictor.Compute2DCoords(
        mol,
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
