import os
from typing import List
from pydantic import BaseModel, Field

from .utils import singleton


@singleton
class ToolSettings(BaseModel):
    ALLOW_LIST: List[str] = Field(default_factory=list, description="List of allowed tools")
    BLOCK_LIST: List[str] = Field(default_factory=list, description="List of blocked tools")
