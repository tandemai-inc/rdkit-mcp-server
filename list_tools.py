import asyncio
from rdkit_mcp.server import mcp
from rdkit_mcp.tools.register_tools import register_tools


async def list_tools():
    """Prints the names and descriptions of all tools registered to the MCP server."""
    print("Registered tools:")
    whitelist = []
    blacklist = []
    await register_tools(mcp, whitelist=whitelist, blacklist=blacklist)
    tool_list = await mcp.list_tools()
    for tool in tool_list:
        name = getattr(tool, 'name', str(tool))
        # desc = getattr(tool, '', '')
        print(f"- {name}")

if __name__ == "__main__":
    asyncio.run(list_tools())
