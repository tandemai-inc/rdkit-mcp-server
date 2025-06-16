from typing import Any, Callable
from mcp.types import (
    AnyFunction,
    ToolAnnotations,
)


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
        fn.tool_enabled = enabled

        # Write module path as Annotation
        tool_module_path = f'{fn.__module__}.{fn.__name__}'
        fn.tool_annotations = {
            "module": ToolAnnotations(title=tool_module_path),
            **(annotations or {}),
        }

        return fn
    return decorator
