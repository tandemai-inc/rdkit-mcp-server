# RDKit MCP Server

🚀 Introducing RDKit Copilot: Agentic Access to RDKit for LLMs
 
🧰 RDKit Copilot is our new open-source MCP server that gives any LLM seamless, agent-level access to RDKit using natural language. No coding required.

🥅 Goal: Expose every function in RDKit 2025.3.1 through an MCP Server
 
🔗 Contribute to making RDKit accessible to LLMs through an open MCP server

💬 Contribute, test, suggest features. We are building it with you.


## Installation

```bash
pip install -e .
```

## Usage

### Start server
```bash
# settings flag optional
python run_server.py --settings=settings.yaml
```

### Start CLI client in another terminal
```bash
export OPENAI_API_KEY="sk-proj-xxx"
python run_client.py
```

## Available Tools
```bash
# settings flag optional
python list_tools.py --settings=settings.yaml
```

