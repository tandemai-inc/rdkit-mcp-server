from ...utils import rdkit_tool, mol_to_sdf, sdf_to_mol

from typing import Tuple, Sequence, Optional, Callable, Any, Iterator, List
from rdkit import Chem
from rdkit.Chem import AllChem, rdchem, rdChemReactions
from rdkit.Geometry import UniformGrid3D


@rdkit_tool()
def assign_bond_orders_from_template(
    ref_sdf_str: str,
    mol_sdf_str: str
) -> str:
    """
    assigns bond orders to a molecule based on the bond orders in a template molecule

    Arguments:
      * refmol: the template molecule
      * mol: the molecule to assign bond orders to
    Note that the template molecule should have no explicit hydrogens else the algorithm will fail.

    It also works if there are different formal charges (this was github issue 235):
    """
    refmol: rdchem.Mol = sdf_to_mol(ref_sdf_str)
    mol: rdchem.Mol = sdf_to_mol(mol_sdf_str)
    result = AllChem.AssignBondOrdersFromTemplate(refmol, mol)
    return mol_to_sdf(result)


@rdkit_tool()
def compute_mol_shape(
    sdf_str: str,
    conf_id: int = -1,
    box_dim: List[int] = None,
    spacing: float = 0.5,
) -> dict:
    """
    returns a grid representation of the molecule's shape

    Parameters:
      * mol (-) - the molecule
      * confId (-) - (optional) the conformer id to use
      * boxDim (-) - (optional) the dimensions of the box to use (Length 3)
      * spacing (-) - (optional) the spacing to use
      * kwargs (-) - additional arguments to pass to the encoding function

    Returns:
    A Dictionary containing information about the UniformGrid3D object
    """
    box_dim = box_dim or [20, 20, 20]
    mol: rdchem.Mol = sdf_to_mol(sdf_str)
    grid: UniformGrid3D = AllChem.ComputeMolShape(mol, conf_id, box_dim, spacing)
    data = grid.GetData()
    dims = (grid.GetXDim(), grid.GetYDim(), grid.GetZDim())
    return {"dims": dims, "spacing": grid.GetSpacing(), "data": data}


@rdkit_tool()
def compute_mol_volume(
    mol_sdf_str: str,
    conf_id: int = -1,
    grid_spacing: float = 0.2,
    box_margin: float = 2.0
) -> float:
    """
    Calculates the volume of a particular conformer of a molecule 

    based on a grid-encoding of the molecular shape.

    Parameters:
      * mol_sdf_str: SDF string of the molecule
      * conf_id (-) - (optional) the conformer id to use
      * grid_spacing (-) - (optional) the spacing to use
      * box_margin (-) - (optional) the margin to use around the molecule
    """
    mol: rdchem.Mol = sdf_to_mol(mol_sdf_str)
    return AllChem.ComputeMolVolume(mol, conf_id, grid_spacing, box_margin)


@rdkit_tool()
def constrained_embed(
    sdf_str: str,
    core_sdf: str,
    use_tethers: bool = True,
    core_conf_id: int = -1,
    random_seed: int = 2342,
    # get_force_field: Optional[Callable[..., Any]] = None,
    **kwargs: Any
) -> Optional[str]:
    """
    generates an embedding with constraints and returns as SDF

    Arguments:
      * sdf_str: SDF string of the molecule to embed
      * core_sdf: SDF string of the template core
      * use_tethers: apply tether forces (optional)
      * core_conf_id: core conformer id (optional)
      * random_seed: RNG seed (optional)
      * get_force_field: custom force field getter (optional)
      * kwargs: additional embedding args
    """
    mol: rdchem.Mol = sdf_to_mol(sdf_str)
    core: rdchem.Mol = sdf_to_mol(core_sdf)
    result = AllChem.ConstrainedEmbed(
        mol,
        core,
        useTethers=use_tethers,
        coreConfId=core_conf_id,
        randomseed=random_seed,
        # getForceField=get_force_field,
        **kwargs
    )
    if result is None:
        return None
    return mol_to_sdf(result)


@rdkit_tool(disabled=True)
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
    mol_sdf_str: str,
    conf_id1: int,
    conf_id2: int,
    atom_ids: Optional[Sequence[int]] = None,
    prealigned: bool = False
) -> float:
    """
    Returns the RMS between two conformations. By default, the conformers will be aligned to the first conformer before the RMS calculation and, as a side-effect, the second will be left in the aligned state.

    Args:
      * sdf_str: SDF string of the molecule
      * conf_id1: first conformer id
      * conf_id2: second conformer id
      * atom_ids: subset of atom indices (optional)
      * prealigned: skip alignment if True (optional)
    """
    mol: rdchem.Mol = sdf_to_mol(mol_sdf_str)
    return AllChem.GetConformerRMS(mol, conf_id1, conf_id2, atomIds=atom_ids, prealigned=prealigned)


@rdkit_tool()
def get_conformer_rms_matrix(
    mol_sdf_str: str,
    atom_ids: Optional[Sequence[int]] = None,
    prealigned: bool = False
) -> List[float]:
    """
    Returns the RMS matrix of the conformers of a molecule. As a side-effect, the conformers will be aligned to the first conformer (i.e. the reference) and will left in the aligned state.

    Parameters:
      * mol_sdf_str (-) - the molecule
      * atomIds (-) - (optional) list of atom ids to use a points for alingment - defaults to all atoms
      * prealigned (-) - (optional) by default the conformers are assumed be unaligned and will therefore be aligned to the first conformer

    Note that the returned RMS matrix is symmetrical, i.e. it is the lower half of the matrix, e.g. for 5 conformers:

        rmsmatrix = [ a,
                      b, c,
                      d, e, f,
                      g, h, i, j]
    where a is the RMS between conformers 0 and 1, b is the RMS between conformers 0 and 2, etc. This way it can be directly used as distance matrix in e.g. Butina clustering.
    """
    mol: rdchem.Mol = sdf_to_mol(mol_sdf_str)
    return AllChem.GetConformerRMSMatrix(mol, atomIds=atom_ids, prealigned=prealigned)


@rdkit_tool()
def transform_mol(
    mol_sdf_str: str,
    transform: Sequence[Sequence[float]],
    conf_id: int = -1,
    keep_confs: bool = False
) -> str:
    """
    Applies the transformation (usually a 4x4 double matrix) to a molecule

    Parameters:
      * mol_sdf_str (-) - the molecule to be transformed
      * tform (-) - the transformation to apply
      * confId (-) - (optional) the conformer id to transform
      * keepConfs (-) - (optional) if keepConfs is False then all but that conformer are removed
    """
    mol: rdchem.Mol = sdf_to_mol(mol_sdf_str)
    result: rdchem.Mol = AllChem.TransformMol(mol, transform, confId=conf_id, keepConfs=keep_confs)
    return mol_to_sdf(result)
