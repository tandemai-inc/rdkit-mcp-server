import os
import sys
import json
import subprocess
import openai
import select
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------
openai.api_key      = os.getenv("OPENAI_API_KEY")
MODEL               = os.getenv("OPENAI_MODEL", "gpt-4o-mini")
MCP_SERVER_CMD      = os.getenv("MCP_SERVER_CMD", "python run_server.py")  # e.g. "python -m my_mcp_server"

# ------------------------------------------------------------------------------
# JSON‑RPC / stdio FRAMEWORK
# ------------------------------------------------------------------------------
# Spawn the MCP server as a subprocess with pipes on stdin/stdout:
process = subprocess.Popen(
    MCP_SERVER_CMD.split(),
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=sys.stderr,  # forward server logs to stderr
    text=True,
    bufsize=0,
)

_rpc_id = 0

def _read_rpc_response():
    # Read headers until blank line
    content_length = None
    while True:

        ready, _, _ = select.select([process.stdout], [], [], 5)
        if not ready:
            logger.warning("No data from MCP server, waiting...")
            continue
        line = process.stdout.readline()
        if not line:
            raise EOFError("MCP server stdout closed")
        line = line.strip()
        if line.lower().startswith("content-length:"):
            content_length = int(line.split(":",1)[1].strip())
        if line == "":
            break
    if content_length is None:
        raise ValueError("No Content-Length header in RPC response")
    # Read the JSON payload
    response = process.stdout.read(content_length)
    return json.loads(response)

def _send_rpc_request(method: str, params: dict):
    global _rpc_id
    _rpc_id += 1

    payload = {
        "jsonrpc": "2.0",
        "id":      _rpc_id,
        "method":  method,
        "params":  params or {}
    }
    body = json.dumps(payload)
    header = f"Content-Length: {len(body.encode('utf-8'))}\r\n\r\n"
    process.stdin.write(header)
    process.stdin.write(body)
    process.stdin.flush()

    resp = _read_rpc_response()
    if "error" in resp:
        raise RuntimeError(f"RPC Error {resp['error']}")
    return resp.get("result")

# ------------------------------------------------------------------------------
# TOOL DISCOVERY & DISPATCH
# ------------------------------------------------------------------------------
def get_exposed_tools():
    """
    RPC call 'list_tools' → expects:
      [ { "name": str,
          "description": str,
          "parameters": <JSON Schema> }, … ]
    """
    return _send_rpc_request("list_tools", {})

def call_tool(tool_name: str, args: dict) -> str:
    """
    RPC call '<tool_name>' with the given args dict,
    returns whatever the server sends back.
    """
    return _send_rpc_request(tool_name, args)

# ------------------------------------------------------------------------------
# MAIN CHAT LOOP
# ------------------------------------------------------------------------------
def main():
    tools = get_exposed_tools()
    functions = [
        {
            "name":        t["name"],
            "description": t.get("description",""),
            "parameters":  t["parameters"],
        }
        for t in tools
    ]

    messages = [
        {"role": "system", "content": "You are a helpful assistant that can call tools when needed."}
    ]

    print(f"🔌 Launched MCP server via stdio with: {MCP_SERVER_CMD}")
    print("Type your queries below (or 'exit').\n")

    while True:
        user_q = input("You: ").strip()
        if user_q.lower() in ("exit","quit"):
            print("👋 Goodbye!")
            break

        messages.append({"role":"user","content":user_q})
        resp = openai.ChatCompletion.create(
            model=MODEL,
            messages=messages,
            functions=functions,
            function_call="auto",
        )
        choice = resp.choices[0].message

        if choice.get("function_call"):
            fname = choice["function_call"]["name"]
            fargs = json.loads(choice["function_call"].get("arguments") or "{}")
            print(f"\n⏳ Calling {fname} with {fargs}…")
            result = call_tool(fname, fargs)
            print(f"✅ {fname} returned:\n{result}\n")

            messages.append({
                "role":"function",
                "name": fname,
                "content": result
            })
            follow = openai.ChatCompletion.create(
                model=MODEL,
                messages=messages
            )
            assistant_reply = follow.choices[0].message["content"]
            print(f"Assistant: {assistant_reply}\n")
            messages.append({"role":"assistant","content":assistant_reply})

        else:
            reply = choice.get("content","")
            print(f"\nAssistant: {reply}\n")
            messages.append({"role":"assistant","content":reply})

if __name__=="__main__":
    main()
