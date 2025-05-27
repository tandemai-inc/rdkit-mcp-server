from ...utils import rdkit_tool

from typing import Tuple, Sequence, Optional, Callable, Any, Iterator, List
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem, rdChemReactions
from rdkit.Geometry import UniformGrid3D


@rdkit_tool()
def assign_bond_orders_from_template(
    refmol: rdchem.Mol,
    mol: rdchem.Mol
) -> rdchem.Mol:
    """
    assigns bond orders to a molecule based on the bond orders in a template molecule

    Arguments:
      * refmol: the template molecule
      * mol: the molecule to assign bond orders to

    Note that the template molecule should have no explicit hydrogens else the algorithm will fail.

    It also works if there are different formal charges (this was github issue 235):
    """
    return AllChem.AssignBondOrdersFromTemplate(refmol, mol)


@rdkit_tool()
def compute_mol_shape(
    mol: rdchem.Mol,
    conf_id: int = -1,
    box_dim: Tuple[int, int, int] = (20, 20, 20),
    spacing: float = 0.5,
    **kwargs: Any
) -> dict:
    """
    returns a grid representation of the molecule’s shape

    Parameters:
      * mol (-) - the molecule
      * confId (-) - (optional) the conformer id to use
      * boxDim (-) - (optional) the dimensions of the box to use
      * spacing (-) - (optional) the spacing to use
      * kwargs (-) - additional arguments to pass to the encoding function

    Returns:
    a UniformGrid3D object
    """
    grid: UniformGrid3D =  AllChem.ComputeMolShape(mol, conf_id, box_dim, spacing, **kwargs)
    data = grid.GetData()          # flat list of floats
    dims = (grid.GetXDim(), grid.GetYDim(), grid.GetZDim())
    return {"dims": dims, "spacing": grid.GetSpacing(), "data": data}


@rdkit_tool()
def compute_mol_volume(
    mol: rdchem.Mol,
    conf_id: int = -1,
    grid_spacing: float = 0.2,
    box_margin: float = 2.0
) -> float:
    """
    Calculates the volume of a particular conformer of a molecule 

    based on a grid-encoding of the molecular shape.

    Parameters:
      * mol (-) - the molecule
      * confId (-) - (optional) the conformer id to use
      * gridSpacing (-) - (optional) the spacing to use
      * boxMargin (-) - (optional) the margin to use around the molecule
    """
    return AllChem.ComputeMolVolume(mol, conf_id, grid_spacing, box_margin)

@rdkit_tool()
def constrained_embed(
    mol: rdchem.Mol,
    core: rdchem.Mol,
    use_tethers: bool = True,
    core_conf_id: int = -1,
    random_seed: int = 2342,
    get_force_field: Optional[Callable[..., Any]] = None,
    **kwargs: Any
) -> Optional[rdchem.Mol]:
    """
    generates an embedding of a molecule where part of the molecule is constrained to have particular coordinates

    Arguments:
      * mol: the molecule to embed
      * core: the molecule to use as a source of constraints
      * useTethers: (optional) if True, the final conformation will be 
optimized subject to a series of extra forces that pull the matching atoms to the positions of the core atoms. Otherwise simple distance constraints based on the core atoms will be used in the optimization.
      * coreConfId: (optional) id of the core conformation to use
      * randomSeed: (optional) seed for the random number generator
      * getForceField: (optional) a function to use to get a force field 
for the final cleanup
      * kwargs: additional arguments to pass to the embedding function
    """
    return AllChem.ConstrainedEmbed(
        mol,
        core,
        useTethers=use_tethers,
        coreConfId=core_conf_id,
        randomseed=random_seed,
        getForceField=get_force_field,
        **kwargs
    )

@rdkit_tool()
def enumerate_library_from_reaction(
    reaction: rdChemReactions.ChemicalReaction,
    sidechain_sets: Sequence[Sequence[rdchem.Mol]],
    return_reactants: bool = False
) -> Iterator[Sequence[rdchem.Mol]]:
    """
    Returns a generator for the virtual library defined by a reaction and a sequence of sidechain sets

    Parameters:
      * reaction (-) - the reaction
      * sidechainSets (-) - a sequence of sequences of sidechains
      * returnReactants (-) - (optional) if True, the generator will return information about the reactants as well as the products
    """
    return AllChem.EnumerateLibraryFromReaction(reaction, sidechain_sets, returnReactants=return_reactants)

@rdkit_tool()
def get_conformer_rms(
    mol: rdchem.Mol,
    conf_id1: int,
    conf_id2: int,
    atom_ids: Optional[Sequence[int]] = None,
    prealigned: bool = False
) -> float:
    """
    Returns the RMS between two conformations. By default, the conformers will be aligned to the first conformer before the RMS calculation and, as a side-effect, the second will be left in the aligned state.

    Parameters:
      * mol (-) - the molecule
      * confId1 (-) - the id of the first conformer
      * confId2 (-) - the id of the second conformer
      * atomIds (-) - (optional) list of atom ids to use a points for alingment - defaults to all atoms
      * prealigned (-) - (optional) by default the conformers are assumed be unaligned and the second conformer be aligned to the first
    """
    return AllChem.GetConformerRMS(mol, conf_id1, conf_id2, atomIds=atom_ids, prealigned=prealigned)

@rdkit_tool()
def get_conformer_rms_matrix(
    mol: rdchem.Mol,
    atom_ids: Optional[Sequence[int]] = None,
    prealigned: bool = False
) -> List[float]:
    """
    Returns the RMS matrix of the conformers of a molecule. As a side-effect, the conformers will be aligned to the first conformer (i.e. the reference) and will left in the aligned state.

    Parameters:
      * mol (-) - the molecule
      * atomIds (-) - (optional) list of atom ids to use a points for alingment - defaults to all atoms
      * prealigned (-) - (optional) by default the conformers are assumed be unaligned and will therefore be aligned to the first conformer

    Note that the returned RMS matrix is symmetrical, i.e. it is the lower half of the matrix, e.g. for 5 conformers:

        rmsmatrix = [ a,
                      b, c,
                      d, e, f,
                      g, h, i, j]
    where a is the RMS between conformers 0 and 1, b is the RMS between conformers 0 and 2, etc. This way it can be directly used as distance matrix in e.g. Butina clustering.
    """
    return AllChem.GetConformerRMSMatrix(mol, atomIds=atom_ids, prealigned=prealigned)

@rdkit_tool()
def transform_mol(
    mol: rdchem.Mol,
    transform: Sequence[Sequence[float]],
    conf_id: int = -1,
    keep_confs: bool = False
) -> rdchem.Mol:
    """
    Applies the transformation (usually a 4x4 double matrix) to a molecule

    Parameters:
      * mol (-) - the molecule to be transformed
      * tform (-) - the transformation to apply
      * confId (-) - (optional) the conformer id to transform
      * keepConfs (-) - (optional) if keepConfs is False then all but that conformer are removed
    """
    return AllChem.TransformMol(mol, transform, confId=conf_id, keepConfs=keep_confs)
