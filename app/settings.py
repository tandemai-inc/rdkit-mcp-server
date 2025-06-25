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


default_file_dir = os.path.join(os.getcwd(), "output")


@singleton
class AppSettings(BaseModel):
    allow_list: List[str] = Field(default_factory=list, description="List of allowed tools")
    block_list: List[str] = Field(default_factory=list, description="List of blocked tools")
    file_dir: DirectoryPath = Field(default=DirectoryPath(default_file_dir), description="Directory where files are stored")
    file_expire: int = Field(default=3600, description="Time in seconds before files expire")


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
