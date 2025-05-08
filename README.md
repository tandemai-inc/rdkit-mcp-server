# RDKit MCP Server

MCP Server providing RDKit cheminformatics tools.

## Installation

```bash
# Installation instructions will go here
# Note: RDKit installation might require specific steps (e.g., via conda)
pip install .
```

## Usage

```bash
# Start server
python run_server.py

# Start client in another terminal window
export OPENAI_API_KEY="sk-proj-xxx"
python run_client.py
```

This will start the server using SSE transport by default.

## Available Tools

- `parse_molecule`: Parses a SMILES string and returns basic properties.
- `draw_molecule`: Generates a PNG image of a molecule from SMILES and returns a file URI.
- `compute_fingerprint`: Computes a molecular fingerprint (default: Morgan).
- `tanimoto_similarity`: Calculates Tanimoto similarity between two molecules.
