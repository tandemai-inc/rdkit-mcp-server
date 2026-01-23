from list_tools import list_tools
import pytest
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


@pytest.mark.asyncio
async def test_list_tools_defaultsettings():
    tool_list = await list_tools()
    assert isinstance(tool_list, list)
    assert len(tool_list) > 0


@pytest.mark.asyncio
async def test_list_tools_allow_list():
    allow_list = [
        "rdkit_mcp.Chem.rdMolDescriptors.CalcPBF"
    ]
    tool_list = await list_tools(allow_list=allow_list)
    assert isinstance(tool_list, list)
    assert len(tool_list) == len(allow_list)


@pytest.mark.asyncio
async def test_list_tools_block_list():
    block_list = [
        "CalcNumUnspecifiedAtomStereoCenters",
        "rdkit_mcp.Chem.rdMolDescriptors.CalcPBF",
        "MolWt"
    ]
    # Get full list of tools
    full_tool_list = await list_tools()
    full_tool_count = len(full_tool_list)

    # Run with block_list, and make sure the count is reduced
    expected_tool_count = full_tool_count - len(block_list)
    tool_list = await list_tools(block_list=block_list)
    assert len(tool_list) == expected_tool_count
