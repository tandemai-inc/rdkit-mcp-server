from typing import Annotated, List
from pydantic import Field
from pydantic import BaseModel


class EncodedFileModel(BaseModel):
    filename: Annotated[str, Field(description="Name of the file.")]
    content: Annotated[str, Field(description="Base 64 encoded bytes containing contents of the file.")]


Smiles = Annotated[str, Field(description="SMILES string representating a molecule's structure")]
MolFragments = Annotated[List[Smiles], Field(description="List of SMILES strings representing molecular fragments")]
PickledMol = Annotated[str, Field(description="Base 64 encoded bytes containing a pickled RDKit Mol object.")]
EncodedFile = Annotated[EncodedFileModel, Field(description="Base 64 encoded bytes containing contents of a file.")]
