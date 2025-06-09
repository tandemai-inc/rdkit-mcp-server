from typing import Annotated
from pydantic import Field


Smiles = Annotated[str, Field(description="SMILES string representating a molecule's structure")]
