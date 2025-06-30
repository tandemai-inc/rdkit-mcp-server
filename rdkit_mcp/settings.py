import logging
import os
from typing import List
from pydantic import BaseModel, Field, DirectoryPath


logger = logging.getLogger(__name__)

default_file_dir = os.path.join(os.getcwd(), "outputs")


def singleton(cls):
    """Add to a class to make it a singleton."""
    instances = {}

    def get_instance(*args, **kwargs):
        if cls not in instances:
            logger.debug("Creating singleton instance of %s with args: %s, kwargs: %s", cls.__name__, args, kwargs)
            instances[cls] = cls(*args, **kwargs)
        elif kwargs and instances[cls] is not None:
            # If kwargs are provided, create a new instance
            # This fixes issue where AppSettings can be instantiated before yaml is loaded,
            # and then we can never get the correct settings.
            logger.debug("Overwriting singleton instance of %s with kwargs: %s", cls.__name__, kwargs)
            instances[cls] = cls(*args, **kwargs)
        return instances[cls]
    return get_instance


@singleton
class AppSettings(BaseModel):
    ALLOW_LIST: List[str] = Field(default_factory=list, description="List of allowed tools")
    BLOCK_LIST: List[str] = Field(default_factory=list, description="List of blocked tools")
    FILE_DIR: DirectoryPath = Field(default=DirectoryPath(default_file_dir), description="Directory where files are stored")
