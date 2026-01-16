"""Custom evaluators for RDKit MCP evaluations."""

import re
from dataclasses import dataclass
from pathlib import Path

from pydantic_ai import Agent, BinaryContent
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


@dataclass
class ImageJudgeEvaluator(Evaluator[TaskInput, TaskOutput]):
    """Evaluator that uses a vision LLM to assess a rendered image.

    Extracts the file path from the task output, reads the image,
    and sends it to a vision-capable model for evaluation.
    """

    rubric: str
    model: str = "openai:gpt-4o"

    async def evaluate(
        self, ctx: EvaluatorContext[TaskInput, TaskOutput]
    ) -> EvaluationReason:
        """Extract image path from output and evaluate with vision LLM."""
        output_text = ctx.output.text

        # Extract file path from output (common patterns)
        path_match = re.search(r"(/[^\s]+\.(?:png|jpg|jpeg|svg|gif))", output_text)
        if not path_match:
            return EvaluationReason(
                value=False,
                reason=f"No image file path found in output: {output_text[:200]}",
            )

        image_path = Path(path_match.group(1))
        if not image_path.exists():
            return EvaluationReason(
                value=False,
                reason=f"Image file not found at path: {image_path}",
            )

        # Determine media type
        suffix = image_path.suffix.lower()
        media_types = {
            ".png": "image/png",
            ".jpg": "image/jpeg",
            ".jpeg": "image/jpeg",
            ".gif": "image/gif",
            ".svg": "image/svg+xml",
        }
        media_type = media_types.get(suffix, "image/png")

        # Read image and create BinaryContent
        image_data = image_path.read_bytes()
        image_content = BinaryContent(data=image_data, media_type=media_type)

        # Build the evaluation prompt
        eval_prompt = f"""You are evaluating a rendered molecule image.

Original task: {ctx.inputs.prompt}

Evaluation rubric: {self.rubric}

Please examine the image and determine if it meets the rubric criteria.
Respond with:
1. PASS or FAIL
2. A brief explanation of your reasoning

Focus on whether the correct substructure is visually highlighted in the image."""

        # Run vision model evaluation
        agent = Agent(self.model)
        result = await agent.run([eval_prompt, image_content])
        response = str(result.output)

        # Parse the response
        passed = "PASS" in response.upper() and "FAIL" not in response.upper()

        return EvaluationReason(
            value=passed,
            reason=response,
        )
