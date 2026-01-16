"""Dataset for comparing batch_map vs synchronous tool execution."""

import csv
from pathlib import Path

from pydantic_evals import Case, Dataset
from pydantic_evals.evaluators import LLMJudge

from evals.evaluators import UsedToolEvaluator
from evals.task import TaskInput

# Common instruction for all cases
REPORT_TOOL_CALLS = "Report all tool calls made with inputs and outputs."


def load_smiles_from_csv(csv_path: str | Path) -> list[str]:
    """Load SMILES strings from the outputs/SMILES.csv file."""
    smiles_list = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            smiles_list.append(row["SMILES"])
    return smiles_list


# Modify DATASET_SIZE to change the number of SMILES used in the dataset
DATASET_SIZE = 500
CSV_PATH = Path(__file__).parent / "data" / "SMILES.csv"
SMILES_FROM_CSV = load_smiles_from_csv(CSV_PATH)[:DATASET_SIZE]

# Build the SMILES list as a comma-separated string for the prompts
SMILES_LIST_STR = ", ".join(SMILES_FROM_CSV)

# Number of SMILES loaded from CSV
SMILES_COUNT = len(SMILES_FROM_CSV)

# Case 1: Synchronous MolWt - call the tool for each SMILES individually
case_sync_molwt = Case(
    name="sync_molwt_big_dataset",
    inputs=TaskInput(
        prompt=(
            f"Calculate the molecular weight for each of these {SMILES_COUNT} SMILES by calling the MolWt tool once for each molecule. "
            f"Do NOT use the batch_map tool - call MolWt individually for each SMILES. "
            f"Response format:\n"
            f"- Don't print smiles and mol weights.\n"
            f"- Confirm the number of molecular weights calculated\n"
            f"Here are the SMILES: {SMILES_LIST_STR}. "
        )
    ),
    expected_output=f"{SMILES_COUNT} molecular weights calculated",
    metadata={
        "smiles_list": SMILES_FROM_CSV,
        "category": "batch_comparison",
        "method": "synchronous",
    },
    evaluators=[
        LLMJudge(
            rubric=(
                f"The output must  state that mol weights for {SMILES_COUNT} molecules were calculated. "
            ),
            include_input=True,
            include_expected_output=True,
        ),
        LLMJudge(
            rubric=(
                f"Confirm that the batch_map tool was NOT used in the response."
            ),
            include_input=True,
            include_expected_output=True,
        ),
    ],
)

# Build the inputs list for batch_map format
BATCH_INPUTS_STR = str([{"smiles": s} for s in SMILES_FROM_CSV])

# Case 2: Batch MolWt using batch_map with high concurrency
case_batch_molwt = Case(
    name="batch_molwt_big_dataset",
    inputs=TaskInput(
        prompt=(
            f"Use the batch_map tool to calculate the molecular weight of {SMILES_COUNT} molecules. "
            f"Use tool_name='MolWt' and pass the inputs as a list of dictionaries with 'smiles' keys. "
            f"Response format:\n"
            f"- Don't print smiles and mol weights.\n"
            f"- Return the number of molecular weights calculated\n"
            f"- Provide a list of each batch_map call, including batch size and concurrency.\n"
            f"Here are the SMILES: {SMILES_LIST_STR}."
        )
    ),
    expected_output=f"{SMILES_COUNT} molecular weights calculated via batch_map with concurrency=100",
    metadata={
        "smiles_list": SMILES_FROM_CSV,
        "category": "batch_comparison",
        "method": "batch_async",
    },
    evaluators=[
        LLMJudge(
            rubric=(
                f"The output must confirm that molecular weights were calculated for all {SMILES_COUNT} molecules. "
                "The batch_map tool must have been used with concurrency set to 100."
            ),
            include_input=True,
            include_expected_output=True,
        ),
        UsedToolEvaluator(tool_name="batch_map"),
    ],
)

# Assemble the dataset
batch_comparison_dataset = Dataset(
    name="batch_comparison_evals",
    cases=[
        case_sync_molwt,
        case_batch_molwt,
    ],
)
