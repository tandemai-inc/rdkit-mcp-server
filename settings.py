import os
from typing import List
from pydantic import BaseModel, Field, DirectoryPath


default_file_dir = os.path.join(os.getcwd(), "outputs")


def singleton(cls):
    """Add to a class to make it a singleton."""
    instances = {}

    def get_instance(*args, **kwargs):
        if cls not in instances:
            instances[cls] = cls(*args, **kwargs)
        return instances[cls]
    return get_instance


@singleton
class AppSettings(BaseModel):
    ALLOW_LIST: List[str] = Field(default_factory=list, description="List of allowed tools")
    BLOCK_LIST: List[str] = Field(default_factory=list, description="List of blocked tools")
    FILE_DIR: DirectoryPath = Field(default=DirectoryPath(default_file_dir), description="Directory where files are stored")
    FILE_EXPIRE: int = Field(default=3600, description="Time in seconds before files expire")


def create_app_settings(yaml_data: dict) -> AppSettings:
    """
    Create an AppSettings instance from a dictionary.

    Args:
        yaml_data (dict): Dictionary containing settings data.

    Returns:
        AppSettings: An instance of AppSettings populated with the provided data.
    """
    return AppSettings(**yaml_data)


def get_app_settings() -> AppSettings:
    """
    Get the singleton instance of AppSettings.

    Returns:
        AppSettings: The singleton instance of AppSettings.
    """
    return AppSettings()
