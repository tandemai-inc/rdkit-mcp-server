from typing import Any, Callable, Optional
from mcp.types import (
    AnyFunction,
    ToolAnnotations,
)
from io import StringIO
from rdkit.Chem import rdchem, SDWriter, ForwardSDMolSupplier


def mol_to_sdf(mol: rdchem.Mol) -> str:
    """
    Helper: Serialize an RDKit Mol object to an SDF string (including properties).

    Args:
      * mol: RDKit molecule to write

    Returns:
      * A complete SDF record string ending with '$$$$'.
    """
    sio = StringIO()
    writer = SDWriter(sio)
    # preserve any arbitrary RDKit properties as SDF tags
    for prop in mol.GetPropNames():
        writer.SetProps({prop: mol.GetProp(prop)})
    writer.write(mol)
    writer.close()
    return sio.getvalue()


def sdf_to_mol(
    sdf_str: str,
    sanitize: bool = True,
    removeHs: bool = True
) -> Optional[rdchem.Mol]:
    """
    Helper: Parse an SDF string into an RDKit Mol object.

    Args:
      * sdf_str: complete SDF record string(s), including '$$$$' delimiters.
      * sanitize: whether to sanitize the molecule upon parsing.
      * removeHs: whether to remove explicit hydrogens.

    Returns:
      * The first RDKit Mol parsed, or None if parsing fails.
    """
    supplier = ForwardSDMolSupplier(StringIO(sdf_str), sanitize=sanitize, removeHs=removeHs)
    mols = list(supplier)
    return mols[0] if mols else None


def rdkit_tool(
    name: str | None = None,
    description: str | None = None,
    annotations: ToolAnnotations | dict[str, Any] | None = None,
    enabled: bool = True,
) -> Callable[[AnyFunction], AnyFunction]:
    """Decorator to register a tool.

    Tools can optionally request a Context object by adding a parameter with the
    Context type annotation. The context provides access to MCP capabilities like
    logging, progress reporting, and resource access.

    Args:
        name: Optional name for the tool (defaults to function name)
        description: Optional description of what the tool does
        annotations: Optional annotations about the tool's behavior
        enabled: If True, the tool will be registered with the MCP server. Useful for tool under development.

    Example:
        @rdkit_tool(name="MyTool", description="This is my tool")
        def my_tool(x: int) -> str:
            return str(x)

        @rdkit_tool(name="MyTool", description="This is my tool")
        def tool_with_context(x: int, ctx: Context) -> str:
            ctx.info(f"Processing {x}")
            return str(x)

        @rdkit_tool(name="MyTool", description="This is my tool")
        async def async_tool(x: int, context: Context) -> str:
            await context.report_progress(50, 100)
            return str(x)
    """

    def decorator(fn: AnyFunction) -> AnyFunction:
        # Add attributes to the function to be used when registering the tool
        fn._is_rdkit_tool = True
        fn.tool_name = name or fn.__name__
        fn.tool_description = description or fn.__doc__
        fn.tool_annotations = annotations or {}
        fn.tool_enabled = enabled
        return fn
    return decorator


def is_rdkit_tool(func: Callable) -> bool:
    """Check if a function is decorated with @rdkit_tool."""
    # Access the original function in case of multiple wrappers
    while hasattr(func, "__wrapped__"):
        func = func.__wrapped__
    return hasattr(func, "_is_rdkit_tool")
