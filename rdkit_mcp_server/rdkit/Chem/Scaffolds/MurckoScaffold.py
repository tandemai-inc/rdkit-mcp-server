from mcp.server.fastmcp.exceptions import ToolError
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from ....tools.utils import rdkit_tool
from ....tools.types import Smiles


@rdkit_tool(description=MurckoScaffold.MakeScaffoldGeneric.__doc__)
def MakeScaffoldGeneric(smiles: Smiles) -> Smiles:
    """Makes a Murcko scaffold generic (i.e. all atom types->C and all bonds ->single."""
    mol = Chem.MolFromSmiles(smiles)
    scaffold: Chem.Mol = MurckoScaffold.MakeScaffoldGeneric(mol)
    if scaffold is None:
        raise ToolError("Failed to generate generic scaffold from the provided SMILES.")
    smiles: Smiles = Chem.MolToSmiles(scaffold)
    return smiles


@rdkit_tool(description=MurckoScaffold.MurckoScaffoldSmilesFromSmiles.__doc__)
def MurckoScaffoldSmilesFromSmiles(smiles: Smiles) -> Smiles:
    """Return the Murcko scaffold as a SMILES string from a SMILES input."""
    return MurckoScaffold.MurckoScaffoldSmilesFromSmiles(smiles)
