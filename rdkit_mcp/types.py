from typing import Annotated
from pydantic import Field


Smiles = Annotated[
    str, Field(description="SMILES string representating a molecule's structure")
]
MolFragments = Annotated[
    tuple[Smiles | None, Smiles],
    Field(
        description="2-tuple (core, side_chain) representing molecular fragments, where core may be a SMILES string of the common scaffold, or `None` if no core could be defined; side_chain is a SMILES string."
    ),
]
PickledMol = Annotated[
    str,
    Field(description="Base 64 encoded bytes containing a pickled RDKit Mol object."),
]
