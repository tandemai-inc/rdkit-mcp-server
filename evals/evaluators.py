"""Custom evaluators for RDKit MCP evaluations."""

from dataclasses import dataclass

from pydantic_evals.evaluators import Evaluator, EvaluatorContext, EvaluationReason

from evals.task import TaskInput, TaskOutput


@dataclass
class UsedToolEvaluator(Evaluator[TaskInput, TaskOutput]):
    """Evaluator that checks if a specific tool was called during task execution."""

    tool_name: str

    def evaluate(
        self, ctx: EvaluatorContext[TaskInput, TaskOutput]
    ) -> EvaluationReason:
        """Check if the required tool was used."""
        tool_names_used = [tc.tool_name for tc in ctx.output.tool_calls]
        tool_was_used = self.tool_name in tool_names_used

        if tool_was_used:
            return EvaluationReason(
                value=True,
                reason=f"Tool '{self.tool_name}' was called. All tools used: {tool_names_used}",
            )
        else:
            return EvaluationReason(
                value=False,
                reason=f"Tool '{self.tool_name}' was NOT called. Tools used: {tool_names_used}",
            )
