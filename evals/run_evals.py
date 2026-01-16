#!/usr/bin/env python
"""CLI entry point for running RDKit MCP evaluations."""

import argparse
import json
import os
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

# Add parent directory to path for direct script execution
if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).parent.parent))

from pydantic_evals import Dataset

from evals.dataset import rdkit_eval_dataset
from evals.task import run_task_sync


def main() -> None:
    # Ensure outputs directory exists
    outputs_dir = Path(os.getcwd()) / "outputs"
    outputs_dir.mkdir(parents=True, exist_ok=True)

    parser = argparse.ArgumentParser(description="Run RDKit MCP Evals")
    parser.add_argument(
        "--output-json",
        type=str,
        help="Path to output JSON results",
    )
    parser.add_argument(
        "--filter",
        type=str,
        help="Filter cases by name (substring match)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed output",
    )
    parser.add_argument(
        "--model",
        type=str,
        default="openai:gpt-4o",
        help="Model to use for evaluation tasks (default: openai:gpt-4o)",
    )
    args = parser.parse_args()

    dataset = rdkit_eval_dataset

    # Apply model override to all cases
    updated_cases = []
    for case in dataset.cases:
        updated_inputs = replace(case.inputs, model=args.model)
        updated_case = replace(case, inputs=updated_inputs)
        updated_cases.append(updated_case)
    dataset = Dataset(name=dataset.name, cases=updated_cases, evaluators=dataset.evaluators)

    # Optional filtering
    if args.filter:
        filtered_cases = [c for c in dataset.cases if args.filter in c.name]
        if not filtered_cases:
            print(f"No cases match filter: {args.filter}")
            sys.exit(1)
        dataset = Dataset(
            name=dataset.name,
            cases=filtered_cases,
        )

    # Run evaluation
    print(f"Running {len(dataset.cases)} evaluation case(s)...")
    report = dataset.evaluate_sync(run_task_sync)

    # Print results
    if args.verbose:
        report.print(include_input=True, include_output=True)
    else:
        report.print()

    # Helper to determine if a case passed (all assertions must have value=True)
    def case_passed(case_result) -> bool:
        if not case_result.assertions:
            return True
        return all(a.value for a in case_result.assertions.values())

    # Output JSON for CI
    if args.output_json:
        case_results: list[dict[str, Any]] = []
        passed_count = 0
        failed_count = 0

        for case_result in report.cases:
            passed = case_passed(case_result)
            if passed:
                passed_count += 1
            else:
                failed_count += 1

            evaluations: list[dict[str, Any]] = []
            for assertion in case_result.assertions.values():
                evaluations.append({
                    "name": assertion.name if hasattr(assertion, "name") else "assertion",
                    "passed": assertion.value if hasattr(assertion, "value") else False,
                    "reason": assertion.reason if hasattr(assertion, "reason") else None,
                })

            case_results.append({
                "name": case_result.name,
                "passed": passed,
                "evaluations": evaluations,
            })

        results = {
            "total": len(report.cases),
            "passed": passed_count,
            "failed": failed_count,
            "cases": case_results,
        }
        with open(args.output_json, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nResults written to {args.output_json}")

    # Exit code for CI
    failed = sum(1 for r in report.cases if not case_passed(r))
    if failed > 0:
        print(f"\n{failed} case(s) failed")
        sys.exit(1)
    else:
        print(f"\nAll {len(report.cases)} case(s) passed")


if __name__ == "__main__":
    main()
