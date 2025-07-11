import os
from typing import List
from pydantic import BaseModel, Field, DirectoryPath


def singleton(cls):
    instances = {}

    def get_instance(*args, **kwargs):
        if cls not in instances:
            instances[cls] = cls(*args, **kwargs)
        return instances[cls]
    return get_instance


default_file_dir = os.path.join(os.getcwd(), "outputs")


@singleton
class AppSettings(BaseModel):
    ALLOW_LIST: List[str] = Field(default_factory=list, description="List of allowed tools")
    BLOCK_LIST: List[str] = Field(default_factory=list, description="List of blocked tools")
    FILE_DIR: DirectoryPath = Field(default=DirectoryPath(default_file_dir), description="Directory where files are stored")
    FILE_EXPIRE: int = Field(default=3600, description="Time in seconds before files expire")
