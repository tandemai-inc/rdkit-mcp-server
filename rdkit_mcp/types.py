from typing import Annotated, List
from pydantic import Field


Smiles = Annotated[str, Field(description="SMILES string representating a molecule's structure")]
MolFragments = Annotated[List[Smiles], Field(description="List of SMILES strings representing molecular fragments")]
