"""Dataset for comparing batch_map vs synchronous tool execution."""

import csv
from pathlib import Path

from pydantic_evals import Case, Dataset
from pydantic_evals.evaluators import LLMJudge

from evals.evaluators import UsedToolEvaluator
from evals.task import TaskInput

# Common instruction for all cases
REPORT_TOOL_CALLS = "Report all tool calls made with inputs and outputs."

DEFAULT_BATCH_SIZE = 100
CSV_PATH = Path(__file__).parent / "data" / "SMILES.csv"


def load_smiles_from_csv(csv_path: str | Path) -> list[str]:
    """Load SMILES strings from the outputs/SMILES.csv file."""
    smiles_list = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            smiles_list.append(row["SMILES"])
    return smiles_list


def create_batch_comparison_dataset(batch_size: int = DEFAULT_BATCH_SIZE) -> Dataset:
    """Create the batch comparison dataset with configurable batch size.

    Args:
        batch_size: Number of SMILES to use from the CSV file.
    """
    smiles_from_csv = load_smiles_from_csv(CSV_PATH)[:batch_size]
    smiles_count = len(smiles_from_csv)

    # Case 1: Synchronous MolWt - call the tool for each SMILES individually
    case_molwt_tool = Case(
        name="sync_molwt_big_dataset",
        inputs=TaskInput(
            prompt=(
                f"Calculate the molecular weight for each of {smiles_count} SMILES by calling the MolWt tool once for each molecule. "
                f"First, call get_smiles_from_context() to retrieve the SMILES list. "
                f"Do NOT use the batch_map tool - call MolWt individually for each SMILES. "
                f"Response format:\n"
                f"- Don't print smiles and mol weights.\n"
                f"- Confirm the number of molecular weights calculated"
            ),
            context={"smiles_list": smiles_from_csv},
        ),
        expected_output=f"{smiles_count} molecular weights calculated",
        metadata={
            "category": "batch_comparison",
            "method": "synchronous",
        },
        evaluators=[
            LLMJudge(
                rubric=(
                    f"The output must state that {smiles_count} mol weights were calculated. "
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

    # Case 2: Batch MolWt using batch_map with high concurrency
    case_batch_map_molwt = Case(
        name="batch_molwt_big_dataset",
        inputs=TaskInput(
            prompt=(
                f"Use the batch_map tool to calculate the molecular weight of {smiles_count} molecules. "
                f"First, call get_smiles_from_context() to retrieve the SMILES list. "
                f"Then pass those SMILES to batch_map with tool_name='MolWt' and inputs as a list of dictionaries with 'smiles' keys. "
                f"Response format:\n"
                f"- Don't print smiles and mol weights.\n"
                f"- Return the number of molecular weights calculated\n"
                f"- Provide a list of each batch_map call, including batch size and concurrency."
            ),
            context={"smiles_list": smiles_from_csv},
        ),
        expected_output=f"{smiles_count} molecular weights calculated via batch_map with concurrency=100",
        metadata={
            "category": "batch_comparison",
            "method": "batch_async",
        },
        evaluators=[
            LLMJudge(
                rubric=(
                    f"The output must confirm that molecular weights were calculated for all {smiles_count} molecules. "
                ),
                include_input=True,
                include_expected_output=True,
            ),
            UsedToolEvaluator(tool_name="batch_map"),
        ],
    )

    # Case 3: compute_descriptors for molecular weight
    case_compute_descriptors_molwt = Case(
        name="compute_descriptors_molwt",
        inputs=TaskInput(
            prompt=(
                f"Use the compute_descriptors tool to calculate the exact molecular weight for {smiles_count} molecules. "
                f"First, call get_smiles_from_context() to retrieve the SMILES list. "
                f"Then call compute_descriptors with descriptor_names=['exactmw']. "
                f"Response format:\n"
                f"- Don't print smiles and mol weights.\n"
                f"- Confirm the number of molecular weights calculated\n"
                f"{REPORT_TOOL_CALLS}"
            ),
            context={"smiles_list": smiles_from_csv},
        ),
        expected_output=f"{smiles_count} molecular weights calculated using compute_descriptors",
        metadata={
            "category": "batch_comparison",
            "method": "compute_descriptors",
        },
        evaluators=[
            LLMJudge(
                rubric=(
                    f"The output must state that {smiles_count} mol weights were calculated. "
                    "The compute_descriptors tool must have been used."
                ),
                include_input=True,
                include_expected_output=True,
            ),
            UsedToolEvaluator(tool_name="compute_descriptors"),
        ],
    )

    return Dataset(
        name="batch_comparison_evals",
        cases=[
            case_molwt_tool,
            case_batch_map_molwt,
            case_compute_descriptors_molwt,
        ],
    )


# Default dataset for backward compatibility
batch_comparison_dataset = create_batch_comparison_dataset()
