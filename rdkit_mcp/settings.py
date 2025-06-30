import os
from typing import List
from pydantic import BaseModel, Field, DirectoryPath

from .utils import singleton


default_file_dir = os.path.join(os.getcwd(), "outputs")


@singleton
class ToolSettings(BaseModel):
    ALLOW_LIST: List[str] = Field(default_factory=list, description="List of allowed tools")
    BLOCK_LIST: List[str] = Field(default_factory=list, description="List of blocked tools")
    FILE_DIR: DirectoryPath = Field(default=DirectoryPath(default_file_dir), description="Directory where files are stored")
