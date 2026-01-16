# RDKit MCP Server: Agentic Access to RDKit for LLMs

RDKit MCP Server is an open-source MCP server that enables language models to interact with RDKit through natural language. The goal is to provide agent-level access to every function in RDKit 2025.3.1 without writing any code.

## Features

* **Seamless Integration**: Exposes RDKit functions via the Model Context Protocol (MCP).
* **Language Model Support**: Connect any LLM that supports the MCP protocol.
* **CLI Client**: Includes a command-line client powered by OpenAI for quick experimentation.

## Table of Contents

* [Installation](#installation)
* [Usage](#usage)

  * [Start the Server](#start-the-server)
  * [CLI Client](#cli-client)
* [Available Tools](#available-tools)
* [Evaluations](#evaluations)
* [Contributing](#contributing)


## Installation

Install the package:

```bash
pip install .
```

## Usage

### Start the Server

```bash
python run_server.py [--settings settings.yaml]
```

See `settings.example.yaml` for setting options

Once the server is running, any MCP-compliant LLM can connect. For example, see the [Claude Desktop quickstart](https://modelcontextprotocol.io/quickstart/user).

### CLI Client

A CLI client is included for rapid prototyping with OpenAI:

```bash
export OPENAI_API_KEY="sk-proj-xxx"
python run_client.py
```

## Available Tools

List all available RDKit tools exposed by the server:

```bash
python list_tools.py [--settings settings.yaml]
```

## Evaluations

The `evals` directory contains a test suite for evaluating RDKit MCP tool outputs and agent responses using [pydantic-evals](https://ai.pydantic.dev/evals/).

### Install Promptfoo

First, install Promptfoo globally using npm:

```bash
pip install ".[evals]"
```

### Start the MCP Server

In one terminal, start the server:

```bash
python run_server.py
```

### Run Evaluations

In another terminal, run the evaluation suite:

```bash
python evals/run_evals.py
```

Options:
- `--verbose` - Show detailed output including inputs and outputs
- `--filter <name>` - Run only cases matching the name
- `--output-json results.json` - Export results to JSON

### Test Cases

The evaluation suite includes 7 test cases:

| Case | Description |
|------|-------------|
| molecular_weight | Calculate molecular weight from SMILES |
| rotatable_bonds | Count rotatable bonds |
| pyridine_detection | Detect pyridine rings in molecules |
| hydrogen_bond_donors_acceptors | Count H-bond donors/acceptors |
| molecular_similarity | Compare molecular similarity |
| atom_highlighting | Render molecule with atom highlighted |
| substructure_highlighting | Render molecule with substructure highlighted |

Each test uses LLM-based evaluation (LLMJudge) to assess whether the agent correctly used the RDKit tools and produced accurate results.

## Contributing

We welcome contributions, feature requests, and bug reports:

See `CONTRIB.md` for guidelines on how to get started.

Together, we can make RDKit accessible to a wider range of applications through natural language interfaces.
