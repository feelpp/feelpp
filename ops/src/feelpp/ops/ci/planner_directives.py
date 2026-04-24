from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import sys
import uuid

SPACK_TARGET = "spack:openmpi"
SPACK_TARGET_ALIASES = {
    "spack": SPACK_TARGET,
    "openmpi": SPACK_TARGET,
    "spack-openmpi": SPACK_TARGET,
}
VALID_SPACK_ONLY_JOBS = {"feelpp-full"}


class PlannerDirectiveError(ValueError):
    """Raised when workflow_dispatch inputs cannot be mapped safely."""


def _split_tokens(raw: str) -> list[str]:
    return [token for token in re.split(r"[\s,]+", raw.strip()) if token]


def _normalize_token_list(raw: str, *, aliases: dict[str, str] | None = None) -> list[str]:
    seen: set[str] = set()
    normalized: list[str] = []
    for token in _split_tokens(raw):
        lower_token = token.lower()
        value = aliases.get(lower_token, lower_token) if aliases else lower_token
        if value in seen:
            continue
        seen.add(value)
        normalized.append(value)
    return normalized


def normalize_targets(raw_targets: str) -> list[str]:
    return _normalize_token_list(raw_targets, aliases=SPACK_TARGET_ALIASES)


def normalize_list_value(raw: str) -> str:
    return ",".join(_normalize_token_list(raw))


def build_planner_message(
    *,
    targets: str = "",
    only: str = "",
    skip: str = "",
    mode: str = "",
) -> str:
    lines: list[str] = []
    normalized_targets = normalize_targets(targets)
    normalized_only = normalize_list_value(only)
    normalized_skip = normalize_list_value(skip)
    normalized_mode = mode.strip().lower()
    contains_spack = SPACK_TARGET in normalized_targets
    spack_only = bool(normalized_targets) and all(target == SPACK_TARGET for target in normalized_targets)

    if normalized_targets:
        lines.append(f"targets={','.join(normalized_targets)}")

    if contains_spack:
        invalid_only = [
            token
            for token in _split_tokens(normalized_only)
            if ":" not in token and token not in VALID_SPACK_ONLY_JOBS
        ]
        if invalid_only:
            invalid_values = ", ".join(invalid_only)
            raise PlannerDirectiveError(
                f"{SPACK_TARGET} only supports the feelpp-full job; invalid only= value(s): {invalid_values}"
            )

        if normalized_skip:
            raise PlannerDirectiveError(
                f"{SPACK_TARGET} does not support skip= filters; use mode=full without component job filters"
            )

        if spack_only and normalized_mode and normalized_mode != "full":
            raise PlannerDirectiveError(f"{SPACK_TARGET} is only supported with mode=full")

        if spack_only and not normalized_mode:
            normalized_mode = "full"

    if normalized_only:
        lines.append(f"only={normalized_only}")
    if normalized_skip:
        lines.append(f"skip={normalized_skip}")
    if normalized_mode:
        lines.append(f"mode={normalized_mode}")

    return "\n".join(lines)


def write_github_output(*, output_path: str, output_name: str, value: str) -> None:
    delimiter = f"FEELPP_OUTPUT_{uuid.uuid4().hex}"
    output_file = Path(output_path)
    with output_file.open("a", encoding="utf-8") as handle:
        handle.write(f"{output_name}<<{delimiter}\n")
        if value:
            handle.write(value)
        handle.write("\n")
        handle.write(f"{delimiter}\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m feelpp.ops.ci.planner_directives")
    parser.add_argument("--targets", default=None, help="Raw workflow_dispatch targets input")
    parser.add_argument("--only", default=None, help="Raw workflow_dispatch only input")
    parser.add_argument("--skip", default=None, help="Raw workflow_dispatch skip input")
    parser.add_argument("--mode", default=None, help="Raw workflow_dispatch mode input")
    parser.add_argument(
        "--github-output",
        default=None,
        help="Path to the GitHub Actions output file. Defaults to $GITHUB_OUTPUT.",
    )
    parser.add_argument(
        "--output-name",
        default="message",
        help="GitHub Actions output variable name to write.",
    )
    return parser


def _resolve_input(value: str | None, env_name: str) -> str:
    if value is not None:
        return value
    return os.environ.get(env_name, "")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        message = build_planner_message(
            targets=_resolve_input(args.targets, "RAW_TARGETS"),
            only=_resolve_input(args.only, "RAW_ONLY"),
            skip=_resolve_input(args.skip, "RAW_SKIP"),
            mode=_resolve_input(args.mode, "RAW_MODE"),
        )
    except PlannerDirectiveError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    output_path = args.github_output or os.environ.get("GITHUB_OUTPUT", "")
    if output_path:
        write_github_output(output_path=output_path, output_name=args.output_name, value=message)
    elif message:
        print(message)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
