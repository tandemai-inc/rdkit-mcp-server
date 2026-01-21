"""Dataset and Case definitions for compute_descriptors tool evaluation."""

from pydantic_evals import Case, Dataset
from pydantic_evals.evaluators import LLMJudge

from evals.evaluators import UsedToolEvaluator
from evals.task import TaskInput

# Common instruction for all cases
REPORT_TOOL_CALLS = "Report all tool calls made with inputs and outputs."


# Case 1: Single molecule, single descriptor
case_single_molecule_single_descriptor = Case(
    name="single_molecule_single_descriptor",
    inputs=TaskInput(
        prompt=(
            "Use the compute_descriptors tool to calculate the exact molecular weight "
            "of ethanol (SMILES: CCO). Use the descriptor name 'exactmw'. "
            f"{REPORT_TOOL_CALLS}"
        )
    ),
    expected_output="~46.042",
    metadata={
        "smiles_list": ["CCO"],
        "descriptor_names": ["exactmw"],
        "category": "basic",
    },
    evaluators=[
        LLMJudge(
            rubric=(
                "The output must report the exact molecular weight of ethanol as approximately 46.04 g/mol. "
                "The compute_descriptors tool must have been used with smiles_list=['CCO'] and "
                "descriptor_names=['exactmw']. Small rounding differences are acceptable."
            ),
            include_input=True,
            include_expected_output=True,
        ),
        UsedToolEvaluator(tool_name="compute_descriptors"),
    ],
)


# Case 2: Single molecule, multiple descriptors
case_single_molecule_multiple_descriptors = Case(
    name="single_molecule_multiple_descriptors",
    inputs=TaskInput(
        prompt=(
            "Use the compute_descriptors tool to calculate the following descriptors for "
            "acetaminophen (SMILES: CC(=O)NC1=CC=C(C=C1)O): exactmw, NumHeavyAtoms, NumHBD, NumHBA. "
            "Report epach descriptor value. "
            f"{REPORT_TOOL_CALLS}"
        )
    ),
    expected_output="exactmw: ~151.06, NumHeavyAtoms: 11, NumHBD: 2, NumHBA: 2",
    metadata={
        "smiles_list": ["CC(=O)NC1=CC=C(C=C1)O"],
        "descriptor_names": ["exactmw", "NumHeavyAtoms", "NumHBD", "NumHBA"],
        "category": "basic",
    },
    evaluators=[
        LLMJudge(
            rubric=(
                "The output must report all four descriptor values for acetaminophen: "
                "1) exactmw approximately 151 g/mol, "
                "2) NumHeavyAtoms = 11, "
                "3) NumHBD (hydrogen bond donors) = 2, "
                "4) NumHBA (hydrogen bond acceptors) = 2. "
                "Small rounding differences for molecular weight are acceptable."
            ),
            include_input=True,
            include_expected_output=True,
        ),
        UsedToolEvaluator(tool_name="compute_descriptors"),
    ],
)


# Case 3: Multiple molecules, single descriptor
case_multiple_molecules_single_descriptor = Case(
    name="multiple_molecules_single_descriptor",
    inputs=TaskInput(
        prompt=(
            "Use the compute_descriptors tool to calculate the number of heavy atoms "
            "(NumHeavyAtoms) for these three molecules in a single call: "
            "ethanol (CCO), benzene (c1ccccc1), and acetic acid (CC(=O)O). "
            "Report the heavy atom count for each molecule. "
            f"{REPORT_TOOL_CALLS}"
        )
    ),
    expected_output="Ethanol: 3, Benzene: 6, Acetic acid: 4",
    metadata={
        "smiles_list": ["CCO", "c1ccccc1", "CC(=O)O"],
        "descriptor_names": ["NumHeavyAtoms"],
        "category": "batch",
    },
    evaluators=[
        LLMJudge(
            rubric=(
                "The output must correctly report NumHeavyAtoms for all three molecules: "
                "ethanol = 3 heavy atoms, benzene = 6 heavy atoms, acetic acid = 4 heavy atoms. "
                "The tool must have processed all molecules in a single call."
            ),
            include_input=True,
            include_expected_output=True,
        ),
        UsedToolEvaluator(tool_name="compute_descriptors"),
    ],
)


# Case 4: Multiple molecules, multiple descriptors (batch computation)
case_batch_multiple_descriptors = Case(
    name="batch_multiple_descriptors",
    inputs=TaskInput(
        prompt=(
            "Use the compute_descriptors tool to compute exactmw and NumAromaticRings "
            "for aspirin (CC(=O)OC1=CC=CC=C1C(=O)O) and caffeine (Cn1cnc2c1c(=O)n(c(=O)n2C)C) "
            "in a single tool call. Report the values for each molecule. "
            f"{REPORT_TOOL_CALLS}"
        )
    ),
    expected_output="Aspirin: MW ~180.04, 1 aromatic ring; Caffeine: MW ~194.08, 0 aromatic rings",
    metadata={
        "smiles_list": ["CC(=O)OC1=CC=CC=C1C(=O)O", "Cn1cnc2c1c(=O)n(c(=O)n2C)C"],
        "descriptor_names": ["exactmw", "NumAromaticRings"],
        "category": "batch",
    },
    evaluators=[
        LLMJudge(
            rubric=(
                "The output must report both descriptors for both molecules: "
                "Aspirin: molecular weight ~180 g/mol and 1 aromatic ring. "
                "Caffeine: molecular weight ~194 g/mol and 0 aromatic rings (the rings in caffeine "
                "are not fully aromatic by RDKit's definition). "
                "Small rounding differences in molecular weight are acceptable."
            ),
            include_input=True,
            include_expected_output=True,
        ),
        UsedToolEvaluator(tool_name="compute_descriptors"),
    ],
)


# Case 5: Drug-likeness screening (Lipinski-related descriptors)
case_drug_likeness_screening = Case(
    name="drug_likeness_screening",
    inputs=TaskInput(
        prompt=(
            "Use the compute_descriptors tool to evaluate drug-likeness properties for "
            "ibuprofen (SMILES: CC(C)CC1=CC=C(C=C1)C(C)C(=O)O). "
            "Calculate these Lipinski-related descriptors: exactmw, CrippenClogP, NumHBD, NumHBA. "
            "Based on Lipinski's Rule of Five (MW < 500, LogP < 5, HBD <= 5, HBA <= 10), "
            "does ibuprofen pass all criteria? "
            f"{REPORT_TOOL_CALLS}"
        )
    ),
    expected_output="MW ~206.13 (<500 ✓), LogP ~3.5 (<5 ✓), HBD 1 (<=5 ✓), HBA 1 (<=10 ✓) - Passes all",
    metadata={
        "smiles_list": ["CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"],
        "descriptor_names": ["exactmw", "CrippenClogP", "NumHBD", "NumHBA"],
        "category": "drug_design",
    },
    evaluators=[
        LLMJudge(
            rubric=(
                "The output must: "
                "1) Report all four descriptor values for ibuprofen, "
                "2) Correctly evaluate each against Lipinski's Rule of Five thresholds, "
                "3) Conclude that ibuprofen passes all Lipinski criteria. "
                "Expected approximate values: MW ~206, LogP ~3-4, HBD 1, HBA ~1-2."
            ),
            include_input=True,
            include_expected_output=True,
        ),
        UsedToolEvaluator(tool_name="compute_descriptors"),
    ],
)


# Case 6: Ring analysis
case_ring_analysis = Case(
    name="ring_analysis",
    inputs=TaskInput(
        prompt=(
            "Use the compute_descriptors tool to analyze the ring structure of naphthalene "
            "(SMILES: c1ccc2ccccc2c1). Calculate NumRings, NumAromaticRings, and NumAliphaticRings. "
            "Explain what each value tells us about the molecule's structure. "
            f"{REPORT_TOOL_CALLS}"
        )
    ),
    expected_output="NumRings: 2, NumAromaticRings: 2, NumAliphaticRings: 0",
    metadata={
        "smiles_list": ["c1ccc2ccccc2c1"],
        "descriptor_names": ["NumRings", "NumAromaticRings", "NumAliphaticRings"],
        "category": "structural",
    },
    evaluators=[
        LLMJudge(
            rubric=(
                "The output must report: "
                "NumRings = 2 (naphthalene has two fused rings), "
                "NumAromaticRings = 2 (both rings are aromatic), "
                "NumAliphaticRings = 0 (no saturated rings). "
                "The explanation should correctly describe naphthalene as having two fused aromatic rings."
            ),
            include_input=True,
            include_expected_output=True,
        ),
        UsedToolEvaluator(tool_name="compute_descriptors"),
    ],
)


# Case 7: Polar surface area comparison
case_tpsa_comparison = Case(
    name="tpsa_comparison",
    inputs=TaskInput(
        prompt=(
            "Use the compute_descriptors tool to compare the topological polar surface area (tpsa) "
            "of three molecules: benzene (c1ccccc1), ethanol (CCO), and glycine (NCC(=O)O). "
            "Which molecule has the highest polar surface area and why? "
            f"{REPORT_TOOL_CALLS}"
        )
    ),
    expected_output="Benzene: 0, Ethanol: ~20, Glycine: ~63. Glycine highest due to amine and carboxylic acid groups.",
    metadata={
        "smiles_list": ["c1ccccc1", "CCO", "NCC(=O)O"],
        "descriptor_names": ["tpsa"],
        "category": "properties",
    },
    evaluators=[
        LLMJudge(
            rubric=(
                "The output must: "
                "1) Report TPSA values for all three molecules (benzene ~0, ethanol ~20, glycine ~60-65), "
                "2) Correctly identify glycine as having the highest TPSA, "
                "3) Explain that glycine's higher TPSA is due to its polar functional groups "
                "(amino and carboxylic acid groups)."
            ),
            include_input=True,
            include_expected_output=True,
        ),
        UsedToolEvaluator(tool_name="compute_descriptors"),
    ],
)


# Case 8: Large batch processing
case_large_batch = Case(
    name="large_batch_processing",
    inputs=TaskInput(
        prompt=(
            "Use the compute_descriptors tool to calculate exactmw and NumRotatableBonds "
            "for these 5 common drugs in a single call: "
            "aspirin (CC(=O)OC1=CC=CC=C1C(=O)O), "
            "acetaminophen (CC(=O)NC1=CC=C(C=C1)O), "
            "ibuprofen (CC(C)CC1=CC=C(C=C1)C(C)C(=O)O), "
            "caffeine (Cn1cnc2c1c(=O)n(c(=O)n2C)C), "
            "and naproxen (CC(C1=CC2=C(C=C1)C=C(C=C2)OC)C(=O)O). "
            "Report the results in a table format. "
            f"{REPORT_TOOL_CALLS}"
        )
    ),
    expected_output="5 drugs with MW and rotatable bonds computed in single call",
    metadata={
        "smiles_list": [
            "CC(=O)OC1=CC=CC=C1C(=O)O",
            "CC(=O)NC1=CC=C(C=C1)O",
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
            "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
            "CC(C1=CC2=C(C=C1)C=C(C=C2)OC)C(=O)O",
        ],
        "descriptor_names": ["exactmw", "NumRotatableBonds"],
        "category": "batch",
    },
    evaluators=[
        LLMJudge(
            rubric=(
                "The output must: "
                "1) Report molecular weight and rotatable bonds for all 5 drugs, "
                "2) Use the compute_descriptors tool with all 5 SMILES in a single call, "
                "3) Present results in a clear tabular or list format. "
                "Expected approximate MWs: aspirin ~180, acetaminophen ~151, ibuprofen ~206, "
                "caffeine ~194, naproxen ~230."
            ),
            include_input=True,
            include_expected_output=True,
        ),
        UsedToolEvaluator(tool_name="compute_descriptors"),
    ],
)


# Assemble the dataset
compute_descriptors_eval_dataset = Dataset(
    name="compute_descriptors_eval_dataset",
    cases=[
        case_single_molecule_single_descriptor,
        case_single_molecule_multiple_descriptors,
        case_multiple_molecules_single_descriptor,
        case_batch_multiple_descriptors,
        case_drug_likeness_screening,
        case_ring_analysis,
        case_tpsa_comparison,
        case_large_batch,
    ],
)
