from __future__ import annotations

import argparse
import json
import sys

from feelpp.ops.common import PREFERRED_DEV_CLI_NAME
from feelpp.pkg.config import discover_repo_root

from .devcontainer import (
    DEFAULT_DEVCONTAINER_TARGET,
    load_profiles,
    validate_generated_files,
    write_generated_files,
)


def _targets_from_args(args: argparse.Namespace) -> list[str] | None:
    if getattr(args, "all", False):
        return None
    return list(args.target) if getattr(args, "target", None) else None


def command_devcontainer_list(args: argparse.Namespace) -> int:
    profiles = load_profiles(args.repo_root, targets=_targets_from_args(args))
    rows = [profile.as_row() for profile in profiles]
    if args.json:
        print(json.dumps(rows, indent=2))
    else:
        for row in rows:
            print(
                f"{row['id']}\t{row['target']}\t{row['cmake_preset']}\t{row['image']}"
            )
    return 0


def _print_result(result) -> None:
    payload = {
        "written": [str(path) for path in result.written],
        "unchanged": [str(path) for path in result.unchanged],
        "drifted": [str(path) for path in result.drifted],
    }
    print(json.dumps(payload, indent=2))


def command_devcontainer_generate(args: argparse.Namespace) -> int:
    result = write_generated_files(
        args.repo_root,
        targets=_targets_from_args(args),
        default_target=args.default_target,
        check=args.check,
    )
    _print_result(result)
    return 1 if result.drifted else 0


def command_devcontainer_validate(args: argparse.Namespace) -> int:
    result = validate_generated_files(
        args.repo_root,
        targets=_targets_from_args(args),
        default_target=args.default_target,
    )
    _print_result(result)
    if result.drifted:
        print("error: generated devcontainer files are out of date", file=sys.stderr)
        return 1
    return 0


def _add_common_generation_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--all",
        action="store_true",
        help="Generate all default Dev Container profiles. This is the default behavior.",
    )
    parser.add_argument(
        "--target",
        action="append",
        default=[],
        help="Restrict generation to an image target from .github/plan-ci.json. Repeatable.",
    )
    parser.add_argument(
        "--default-target",
        default=DEFAULT_DEVCONTAINER_TARGET,
        help="Image target used for .devcontainer/devcontainer.json.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=PREFERRED_DEV_CLI_NAME)
    parser.add_argument(
        "--repo-root",
        default=str(discover_repo_root()),
        help="Path to the Feel++ repository root",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    devcontainer_parser = subparsers.add_parser(
        "devcontainer",
        help="Generate and validate Dev Container profiles",
    )
    devcontainer_subparsers = devcontainer_parser.add_subparsers(
        dest="devcontainer_command",
        required=True,
    )

    list_parser = devcontainer_subparsers.add_parser(
        "list",
        help="List generated Dev Container profiles",
    )
    list_parser.add_argument(
        "--target",
        action="append",
        default=[],
        help="Restrict output to an image target from .github/plan-ci.json. Repeatable.",
    )
    list_parser.add_argument("--json", action="store_true", help="Print JSON output")
    list_parser.set_defaults(func=command_devcontainer_list)

    generate_parser = devcontainer_subparsers.add_parser(
        "generate",
        help="Write generated Dev Container profiles and CMake user presets",
    )
    _add_common_generation_args(generate_parser)
    generate_parser.add_argument(
        "--check",
        action="store_true",
        help="Report drift without writing files",
    )
    generate_parser.set_defaults(func=command_devcontainer_generate)

    validate_parser = devcontainer_subparsers.add_parser(
        "validate",
        help="Check that generated Dev Container profiles are up to date",
    )
    _add_common_generation_args(validate_parser)
    validate_parser.set_defaults(func=command_devcontainer_validate)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
