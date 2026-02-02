import logging

from pydantic import BaseModel, Field
from rdkit.Chem import Mol
from rdkit.Chem import rdDistGeom as _rdDistGeom
from ..decorators import rdkit_tool
from ..types import PickledMol
from ..utils import decode_mol, encode_mol

logger = logging.getLogger(__name__)


class EmbedParameters(BaseModel):
    """Parameters for controlling molecule embedding with distance geometry.

    This schema mirrors the RDKit EmbedParameters class and provides configuration
    options for 3D coordinate generation.
    """
    # Core Embedding Control
    maxIterations: int = Field(
        default=0,
        description="Maximum number of embedding attempts per conformation. 0 means use default (10*numAtoms)"
    )
    randomSeed: int = Field(
        default=-1,
        description="Seed for random number generator. Use -1 for non-deterministic behavior"
    )
    clearConfs: bool = Field(
        default=True,
        description="Remove existing conformations before embedding"
    )
    useRandomCoords: bool = Field(
        default=False,
        description="Start from random coordinates instead of eigenvalues"
    )
    boxSizeMult: float = Field(
        default=2.0,
        description="Controls random coordinate box size. Positive: multiplies largest distance; Negative: absolute size"
    )

    # Eigenvalue Handling
    randNegEig: bool = Field(
        default=True,
        description="Pick coordinates at random if negative eigenvalues occur"
    )
    numZeroFail: int = Field(
        default=1,
        description="Fail embedding if zero eigenvalues reach this threshold"
    )

    # Minimization & Optimization
    optimizerForceTol: float = Field(
        default=0.001,
        description="Tolerance for distance-geometry force field minimization"
    )
    basinThresh: float = Field(
        default=5.0,
        description="Basin threshold for distance geometry force field"
    )
    boundsMatForceScaling: float = Field(
        default=1.0,
        description="Weight scaling for atom pair distance restraints"
    )

    # Structural Constraints
    enforceChirality: bool = Field(
        default=True,
        description="Enforce correct chirality at stereogenic centers"
    )
    useExpTorsionAnglePrefs: bool = Field(
        default=True,
        description="Apply experimental torsion angle preferences"
    )
    useBasicKnowledge: bool = Field(
        default=True,
        description="Impose basic-knowledge constraints such as flat rings"
    )
    useSmallRingTorsions: bool = Field(
        default=False,
        description="Apply small ring torsion preferences"
    )
    useMacrocycleTorsions: bool = Field(
        default=True,
        description="Apply macrocycle-specific torsion angles"
    )
    useMacrocycle14config: bool = Field(
        default=True,
        description="Use ETKDGv3 1-4 distance bounds for macrocycles"
    )
    forceTransAmides: bool = Field(
        default=True,
        description="Constrain amide bonds to trans configuration"
    )

    # Multi-Conformer Generation
    pruneRmsThresh: float = Field(
        default=-1.0,
        description="Keep conformations at least this distance apart (RMSD). Use -1.0 to disable pruning"
    )
    onlyHeavyAtomsForRMS: bool = Field(
        default=False,
        description="Use only heavy atoms for RMSD filtering"
    )
    useSymmetryForPruning: bool = Field(
        default=True,
        description="Apply molecular symmetry during pruning"
    )
    symmetrizeConjugatedTerminalGroupsForPruning: bool = Field(
        default=True,
        description="Symmetrize terminal conjugated groups for pruning"
    )
    numThreads: int = Field(
        default=1,
        description="Number of threads for multiple conformation generation. 0 uses all available cores"
    )

    # Advanced Options
    ETversion: int = Field(
        default=2,
        description="Experimental torsion definition version (1 or 2). Use 2 for both ETKDGv2 and ETKDGv3"
    )
    ignoreSmoothingFailures: bool = Field(
        default=False,
        description="Continue embedding despite bounds matrix smoothing failure"
    )
    timeout: int = Field(
        default=0,
        description="Maximum seconds per fragment. 0 means unlimited"
    )
    trackFailures: bool = Field(
        default=False,
        description="Monitor failure points during embedding"
    )
    embedFragmentsSeparately: bool = Field(
        default=True,
        description="Process disconnected fragments independently"
    )
    enableSequentialRandomSeeds: bool = Field(
        default=False,
        description="Support conformer generation restart with sequential seeds"
    )
    verbose: bool = Field(
        default=False,
        description="Enable detailed output during embedding"
    )

    # Coordinate Map (for fixed atom positions)
    coordMap: dict[int, tuple[float, float, float]] = Field(
        default_factory=dict,
        description="Dictionary mapping atom IDs to fixed coordinates (x, y, z)"
    )

    def to_rdkit_params(self, method: str = "ETKDGv3"):
        """Convert this Pydantic model to an RDKit EmbedParameters object.

        Args:
            method: The RDKit parameter initialization method to use.
                    Currently supports: "ETKDGv3" (default)
                    Future support could include: "ETKDGv2", "ETKDG", "KDG", "srETKDGv3", etc.

        Returns:
            Configured RDKit EmbedParameters object
        """
        # Create the appropriate RDKit params object
        if method == "ETKDGv3":
            rdkit_params = _rdDistGeom.ETKDGv3()
        else:
            raise ValueError(f"Unsupported method: {method}. Currently only 'ETKDGv3' is supported.")

        # Set all parameters from this Pydantic model
        params_dict = self.model_dump()
        coord_map = params_dict.pop('coordMap', None)

        for param_name, param_value in params_dict.items():
            if hasattr(rdkit_params, param_name):
                setattr(rdkit_params, param_name, param_value)
            else:
                logger.warning(f"Parameter '{param_name}' not found on RDKit EmbedParameters object")

        # Set coordMap if provided (requires special method)
        if coord_map:
            rdkit_params.SetCoordMap(coord_map)

        return rdkit_params


class EmbedMoleculeResult(BaseModel):
    """Result from embedding a molecule with 3D coordinates.

    Attributes:
        conf_id: ID of the new conformation added to the molecule.
                 Returns -1 if the embedding fails.
        mol: The modified molecule with embedded 3D coordinates (pickled).
    """
    conf_id: int = Field(
        description="ID of the new conformation added to the molecule, or -1 if embedding fails"
    )
    mol: PickledMol = Field(
        description="The molecule with embedded 3D coordinates"
    )


class EmbedMultipleConfsResult(BaseModel):
    """Result from embedding multiple conformations of a molecule.

    Attributes:
        conf_ids: List of IDs for the new conformations added to the molecule.
                  Empty list if embedding fails completely.
        num_confs: Number of conformations successfully embedded.
        mol: The modified molecule with embedded 3D coordinates (pickled).
    """
    conf_ids: list[int] = Field(
        description="List of conformation IDs added to the molecule"
    )
    num_confs: int = Field(
        description="Number of conformations successfully embedded"
    )
    mol: PickledMol = Field(
        description="The molecule with embedded 3D coordinates for all conformations"
    )


@rdkit_tool(description=_rdDistGeom.EmbedMolecule.__doc__)
def EmbedMolecule(
    p_mol: PickledMol,
    params: EmbedParameters = None
) -> EmbedMoleculeResult:
    """
    Use distance geometry to obtain initial coordinates for a molecule.

    ARGUMENTS:
        p_mol: the pickled molecule of interest
        params: EmbedParameters object controlling the embedding behavior.
            If not provided, default parameters will be used.

    RETURNS:
        EmbedMoleculeResult with:
            - conf_id: ID of the new conformation added to the molecule or -1 if the embedding fails
            - mol: The modified molecule with embedded 3D coordinates (pickled)
    """
    mol: Mol = decode_mol(p_mol)

    # Use default parameters if none provided
    if params is None:
        params = EmbedParameters()

    # Convert Pydantic parameters to RDKit parameters
    rdkit_params: _rdDistGeom.EmbedParameters = params.to_rdkit_params()

    conf_id = _rdDistGeom.EmbedMolecule(mol, rdkit_params)

    return EmbedMoleculeResult(
        conf_id=conf_id,
        mol=encode_mol(mol)
    )


@rdkit_tool(description=_rdDistGeom.EmbedMultipleConfs.__doc__)
def EmbedMultipleConfs(
    p_mol: PickledMol,
    numConfs: int = 10,
    params: EmbedParameters = None
) -> EmbedMultipleConfsResult:
    """
    Use distance geometry to obtain multiple sets of coordinates for a molecule.

    ARGUMENTS:
        p_mol: the pickled molecule of interest
        numConfs: the number of conformers to generate (default: 10)
        params: EmbedParameters object controlling the embedding behavior.
            If not provided, default parameters will be used.
            Note: pruneRmsThresh in params controls conformer diversity - only conformations
            at least this distance apart (RMSD in Angstroms) will be retained.

    RETURNS:
        EmbedMultipleConfsResult with:
            - conf_ids: List of conformation IDs added to the molecule
            - num_confs: Number of conformations successfully embedded
            - mol: The modified molecule with all embedded conformations (pickled)
    """
    mol: Mol = decode_mol(p_mol)

    # Use default parameters if none provided
    if params is None:
        params = EmbedParameters()

    # Convert Pydantic parameters to RDKit parameters
    rdkit_params: _rdDistGeom.EmbedParameters = params.to_rdkit_params()

    # Embed multiple conformations
    conf_ids = _rdDistGeom.EmbedMultipleConfs(mol, numConfs, rdkit_params)

    return EmbedMultipleConfsResult(
        conf_ids=list(conf_ids),
        num_confs=len(conf_ids),
        mol=encode_mol(mol)
    )
