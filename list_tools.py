import asyncio
from register_tools import register_tools
from run_server import mcp


async def list_tools():
    """Prints the names and descriptions of all tools registered to the MCP server."""
    print("Registered tools:")
    whitelist = []
    blacklist = []
    await register_tools(mcp, whitelist=whitelist, blacklist=blacklist)
    tool_list = await mcp.list_tools()
    for tool in tool_list:
        module_path = tool.annotations.module.title
        print(f"- {module_path}")

if __name__ == "__main__":
    asyncio.run(list_tools())
