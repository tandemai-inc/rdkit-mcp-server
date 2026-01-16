"""Dataset and Case definitions for RDKit MCP evaluations."""

from pydantic_evals import Case, Dataset
from pydantic_evals.evaluators import LLMJudge

from evals.evaluators import ImageJudgeEvaluator, UsedToolEvaluator
from evals.task import TaskInput

# Case 1: Molecular Weight
case_molecular_weight = Case(
    name="molecular_weight",
    inputs=TaskInput(
        prompt="What is the molecular weight of the compound with SMILES: CC(=O)NC1=CC=C(C=C1)O?"
    ),
    expected_output="151.165 g/mol",
    metadata={"smiles": "CC(=O)NC1=CC=C(C=C1)O", "category": "properties"},
    evaluators=[
        LLMJudge(
            rubric="The final output states the Molecular weight is 151.165 g/mol. Rounding is acceptable.",
            include_input=True,
            include_expected_output=True,
        )
    ],
)

# Case 2: Rotatable Bonds
case_rotatable_bonds = Case(
    name="rotatable_bonds",
    inputs=TaskInput(
        prompt="Give me the number of rotatable bonds in this SMILES: CCOC(=O)C1=CC=CC=C1?"
    ),
    expected_output="2",
    metadata={"smiles": "CCOC(=O)C1=CC=CC=C1", "category": "properties"},
    evaluators=[
        LLMJudge(
            rubric="The final output states the number of rotatable bonds is 2.",
            include_input=True,
        ),
        LLMJudge(
            rubric="The response should contain an output type which returns the value 2.",
            include_input=True,
        ),
    ],
)

# Case 3: Pyridine Ring Detection
case_pyridine_detection = Case(
    name="pyridine_detection",
    inputs=TaskInput(
        prompt="Find all molecules in this list that contain a pyridine ring: C1=CC=CN=C1, C1=CC=CC=C1, CC(C)C1=CN=CC=C1"
    ),
    expected_output="C1=CC=CN=C1, CC(C)C1=CN=CC=C1",
    metadata={"category": "substructure"},
    evaluators=[
        LLMJudge(
            rubric="The final output states that C1=CC=CN=C1 and CC(C)C1=CN=CC=C1 contain pyridine rings.",
            include_input=True,
            include_expected_output=True,
        )
    ],
)

# Case 4: Hydrogen Bond Donors/Acceptors
case_hbond = Case(
    name="hydrogen_bond_donors_acceptors",
    inputs=TaskInput(
        prompt="Compute the number of hydrogen bond donors and acceptors for this molecule: COc1ccc(CO)cc1"
    ),
    expected_output="Donors: 1, Acceptors: 2",
    metadata={"smiles": "COc1ccc(CO)cc1", "category": "properties"},
    evaluators=[
        LLMJudge(
            rubric="The final output states that the SMILES has 1 Donor hydrogen bond, and 2 Acceptors.",
            include_input=True,
        )
    ],
)

# Case 5: Molecular Similarity
case_similarity = Case(
    name="molecular_similarity",
    inputs=TaskInput(
        prompt="""Which of these is more similar to this reference molecule
Reference: CC(=O)OC1=CC=CC=C1C(=O)O
Candidates: CC(OC1C(C(O)=O)=CC(C)=CC=1)=O or CCOC(C)C or C1C=CC=CC=1"""
    ),
    expected_output="CC(OC1C(C(O)=O)=CC(C)=CC=1)=O",
    metadata={"category": "similarity"},
    evaluators=[
        LLMJudge(
            rubric="The final output states that CC(OC1C(C(O)=O)=CC(C)=CC=1)=O is the most similar to the reference",
            include_input=True,
        )
    ],
)

# Case 6: Atom Highlighting
case_atom_highlight = Case(
    name="atom_highlighting",
    inputs=TaskInput(
        prompt="Render an image of this SMILES, and highlight atom 1: CC(=O)NC1=CC=CC=C1"
    ),
    expected_output="Image with atom 1 highlighted",
    metadata={"smiles": "CC(=O)NC1=CC=CC=C1", "category": "rendering"},
    evaluators=[
        LLMJudge(
            rubric="The final output should render an image with atom 1 highlighted in the molecule CC(=O)NC1=CC=CC=C1.",
            include_input=True,
        )
    ],
)

# Case 7: Substructure Highlighting
case_substructure_highlight = Case(
    name="substructure_highlighting",
    inputs=TaskInput(
        prompt="Highlight the amide substructure in this compound: CC(=O)NC1=CC=CC=C1"
    ),
    expected_output="Image with amide highlighted",
    metadata={"smiles": "CC(=O)NC1=CC=CC=C1", "category": "rendering"},
    evaluators=[
        LLMJudge(
            rubric="The final output should render an image with the amide substructure highlighted in the molecule CC(=O)NC1=CC=CC=C1.",
            include_input=True,
        ),
        ImageJudgeEvaluator(
            rubric="The image shows a molecule with the amide functional group (C(=O)N, carbonyl bonded to nitrogen) clearly highlighted in a distinct color. The highlighting should visually distinguish the amide substructure from the rest of the molecule.",
        ),
    ],
)

# Case 8: Batch Molecular Weight Calculation
case_batch_molecular_weight = Case(
    name="batch_molecular_weight",
    inputs=TaskInput(
        prompt=(
            "Use the batch_map tool to calculate the molecular weight of the following 3 molecules in a single batch call: "
            "CCO (ethanol), CC(=O)O (acetic acid), and C1=CC=CC=C1 (benzene). "
            "Report the molecular weight for each molecule."
            "Report all tool calls made with inputs and outputs"
            "I expect the batch_map tool to be called"
        )
    ),
    expected_output="Ethanol: 46.069, Acetic acid: 60.052, Benzene: 78.114",
    metadata={
        "smiles_list": ["CCO", "CC(=O)O", "C1=CC=CC=C1"],
        "category": "batch",
    },
    evaluators=[
        LLMJudge(
            rubric=(
                "The output must report molecular weights for all 3 molecules. "
                "Expected approximate values: ethanol ~46 g/mol, acetic acid ~60 g/mol, benzene ~78 g/mol. "
                "Small rounding differences are acceptable."
            ),
            include_input=True,
            include_expected_output=True,
        ),
        UsedToolEvaluator(tool_name="batch_map"),
    ],
)

# Assemble the dataset
rdkit_eval_dataset = Dataset(
    name="rdkit_mcp_evals",
    cases=[
        case_molecular_weight,
        case_rotatable_bonds,
        case_pyridine_detection,
        case_hbond,
        case_similarity,
        case_atom_highlight,
        case_substructure_highlight,
        case_batch_molecular_weight,
    ],
    
)
