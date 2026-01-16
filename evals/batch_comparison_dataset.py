"""Dataset for comparing batch_map vs synchronous tool execution."""

from pydantic_evals import Case, Dataset
from pydantic_evals.evaluators import LLMJudge

from evals.evaluators import UsedToolEvaluator
from evals.task import TaskInput

# Common instruction for all cases
REPORT_TOOL_CALLS = "Report all tool calls made with inputs and outputs."

# 100 diverse SMILES for batch comparison tests
SMILES_100 = [
    "CCO",                      # ethanol
    "CC(=O)O",                  # acetic acid
    "C1=CC=CC=C1",              # benzene
    "CCN(CC)CC",                # triethylamine
    "CC(C)O",                   # isopropanol
    "CCCC",                     # butane
    "CC(=O)C",                  # acetone
    "CCCCCC",                   # hexane
    "c1ccccc1O",                # phenol
    "CC(C)(C)O",                # tert-butanol
    "CCOCC",                    # diethyl ether
    "CCN",                      # ethylamine
    "CCC(=O)O",                 # propionic acid
    "CCCO",                     # propanol
    "CC=C",                     # propene
    "C#N",                      # hydrogen cyanide
    "CC#N",                     # acetonitrile
    "CCCCCCCC",                 # octane
    "c1ccc(O)cc1O",             # catechol
    "CC(=O)OCC",                # ethyl acetate
    "CCCCO",                    # butanol
    "CCCCCO",                   # pentanol
    "c1ccc(N)cc1",              # aniline
    "CC(C)CC",                  # isopentane
    "CCC(C)C",                  # isopentane isomer
    "CCOC(=O)C",                # ethyl acetate isomer
    "c1ccc(Cl)cc1",             # chlorobenzene
    "c1ccc(Br)cc1",             # bromobenzene
    "c1ccc(F)cc1",              # fluorobenzene
    "c1ccc(I)cc1",              # iodobenzene
    "CC(=O)N",                  # acetamide
    "CCCCCCCCCC",               # decane
    "c1cccnc1",                 # pyridine
    "c1ccncc1",                 # pyridine isomer
    "CC(C)C(=O)O",              # isobutyric acid
    "CCCCC(=O)O",               # valeric acid
    "c1ccc2ccccc2c1",           # naphthalene
    "CCc1ccccc1",               # ethylbenzene
    "Cc1ccccc1",                # toluene
    "Cc1ccc(C)cc1",             # p-xylene
    "Cc1cccc(C)c1",             # m-xylene
    "Cc1ccccc1C",               # o-xylene
    "c1ccc(C(C)C)cc1",          # cumene
    "CCCCCCCCCCC",              # undecane
    "CCCCCCCCCCCC",             # dodecane
    "CC(=O)Nc1ccccc1",          # acetanilide
    "OCC(O)CO",                 # glycerol
    "OCCO",                     # ethylene glycol
    "OCCCO",                    # propylene glycol
    "CC(O)C",                   # isopropanol duplicate check
    "CCC(C)(C)O",               # tert-amyl alcohol
    "c1ccc(OC)cc1",             # anisole
    "COc1ccc(O)cc1",            # 4-methoxyphenol
    "c1ccc(C=O)cc1",            # benzaldehyde
    "CC(=O)c1ccccc1",           # acetophenone
    "c1ccc(C(=O)O)cc1",         # benzoic acid
    "c1ccc(C(=O)N)cc1",         # benzamide
    "CCCCCl",                   # 1-chlorobutane
    "CCCCBr",                   # 1-bromobutane
    "CCCCI",                    # 1-iodobutane
    "CCCCF",                    # 1-fluorobutane
    "CC(Cl)C",                  # 2-chloropropane
    "CC(Br)C",                  # 2-bromopropane
    "CCC(=O)CC",                # 3-pentanone
    "CCCC(=O)C",                # 2-pentanone
    "CC(=O)CCC",                # 2-pentanone isomer
    "c1ccc(N(C)C)cc1",          # N,N-dimethylaniline
    "c1ccc(NC)cc1",             # N-methylaniline
    "NCCN",                     # ethylenediamine
    "NCCCN",                    # 1,3-diaminopropane
    "NCCCCN",                   # 1,4-diaminobutane
    "OC(=O)C(O)=O",             # oxalic acid
    "OC(=O)CC(=O)O",            # malonic acid
    "OC(=O)CCC(=O)O",           # succinic acid
    "OC(=O)CCCC(=O)O",          # glutaric acid
    "OC(=O)CCCCC(=O)O",         # adipic acid
    "c1ccc2ncccc2c1",           # quinoline
    "c1cnc2ccccc2n1",           # quinazoline
    "c1ccc2[nH]ccc2c1",         # indole
    "c1ccc2occc2c1",            # benzofuran
    "c1ccc2sccc2c1",            # benzothiophene
    "C1CCCCC1",                 # cyclohexane
    "C1CCCC1",                  # cyclopentane
    "C1CCC1",                   # cyclobutane
    "C1CC1",                    # cyclopropane
    "C1CCCCCCC1",               # cyclooctane
    "c1ccc(CCO)cc1",            # 2-phenylethanol
    "c1ccc(CCCO)cc1",           # 3-phenylpropanol
    "CCCCCCC(=O)O",             # heptanoic acid
    "CCCCCCCC(=O)O",            # octanoic acid
    "CCCCCCCCC(=O)O",           # nonanoic acid
    "CCCCCCCCCC(=O)O",          # decanoic acid
    "c1ccc(S)cc1",              # thiophenol
    "CCSC",                     # ethyl methyl sulfide
    "CCSCC",                    # diethyl sulfide
    "CS(=O)C",                  # dimethyl sulfoxide
    "CC(=O)SC",                 # methyl thioacetate
    "c1ccoc1",                  # furan
    "c1ccsc1",                  # thiophene
    "c1cc[nH]c1",               # pyrrole
]

# Build the SMILES list as a comma-separated string for the prompts
SMILES_LIST_STR = ", ".join(SMILES_100)

# Case 1: Synchronous MolWt - call the tool 100 times individually
case_sync_molwt_100 = Case(
    name="sync_molwt_100",
    inputs=TaskInput(
        prompt=(
            f"Calculate the molecular weight for each of these 100 SMILES by calling the MolWt tool once for each molecule. "
            f"Do NOT use the batch_map tool - call MolWt individually for each SMILES. "
            f"Here are the SMILES: {SMILES_LIST_STR}. "
            f"Report the count of molecules processed and confirm all 100 were calculated. "
            f"{REPORT_TOOL_CALLS}"
        )
    ),
    expected_output="100 molecular weights calculated",
    metadata={
        "smiles_list": SMILES_100,
        "category": "batch_comparison",
        "method": "synchronous",
    },
    evaluators=[
        LLMJudge(
            rubric=(
                "The output must confirm that molecular weights were calculated for all 100 molecules. "
                "The response should indicate that MolWt was called individually (not via batch_map)."
            ),
            include_input=True,
            include_expected_output=True,
        ),
    ],
)

# Build the inputs list for batch_map format
BATCH_INPUTS_STR = str([{"smiles": s} for s in SMILES_100])

# Case 2: Batch MolWt using batch_map with concurrency=100
case_batch_molwt_100 = Case(
    name="batch_molwt_100",
    inputs=TaskInput(
        prompt=(
            f"Use the batch_map tool to calculate the molecular weight of 100 molecules in a single batch call. "
            f"Set the concurrency parameter to 100 to run all calculations in parallel. "
            f"Use tool_name='MolWt' and pass the inputs as a list of dictionaries with 'smiles' keys. "
            f"Here are the SMILES: {SMILES_LIST_STR}. "
            f"Report the count of molecules processed and confirm all 100 were calculated. "
            f"{REPORT_TOOL_CALLS}"
        )
    ),
    expected_output="100 molecular weights calculated via batch_map with concurrency=100",
    metadata={
        "smiles_list": SMILES_100,
        "category": "batch_comparison",
        "method": "batch_async",
        "concurrency": 100,
    },
    evaluators=[
        LLMJudge(
            rubric=(
                "The output must confirm that molecular weights were calculated for all 100 molecules. "
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
        case_sync_molwt_100,
        case_batch_molwt_100,
    ],
)
