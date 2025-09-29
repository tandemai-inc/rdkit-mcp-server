import os

from fastapi import Request
from fastapi.responses import FileResponse, JSONResponse
from mcp.server.fastmcp import FastMCP

from rdkit_mcp.settings import ToolSettings

__all__ = ["app"]


app = FastMCP("RDKit-MCP Server")


@app.custom_route("/files", methods=["GET"])
async def list_files(request: Request):
    settings = ToolSettings()
    file_dir = str(settings.FILE_DIR)
    try:
        files = os.listdir(file_dir)
        return JSONResponse({"files": files})
    except Exception as e:
        return JSONResponse({"error": str(e)})


@app.custom_route("/files/{filename}", methods=["GET"])
async def get_file(request: Request):
    filename = request.path_params.get("filename")
    settings = ToolSettings()
    file_dir = str(settings.FILE_DIR)
    file_path = os.path.join(file_dir, filename)
    if not os.path.isfile(file_path):
        return JSONResponse(status_code=404, content={"error": "File not found"})
    return FileResponse(path=file_path, filename=filename)
