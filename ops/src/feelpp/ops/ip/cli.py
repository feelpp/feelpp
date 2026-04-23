from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from feelpp.ops.common import PREFERRED_IP_CLI_NAME
from feelpp.pkg.config import discover_repo_root

from .exporter import build_public_bundle, write_public_bundle
from .metrics import collect_code_metrics
from .repository import PublicMetadataRepository
from .yamlio import dump_yaml


def _repository_from_args(args: argparse.Namespace) -> PublicMetadataRepository:
    return PublicMetadataRepository(
        repo_root=args.repo_root,
        metadata_path=args.metadata_path,
    )


def command_validate_public(args: argparse.Namespace) -> int:
    repository = _repository_from_args(args)
    repository.read_public_metadata()
    print(f"{repository.relative_path(repository.metadata_path)}: ok")
    return 0


def command_show_public(args: argparse.Namespace) -> int:
    repository = _repository_from_args(args)
    payload = repository.read_public_metadata()
    if args.format == "json":
        print(json.dumps(payload, indent=2, sort_keys=False))
    else:
        print(dump_yaml(payload), end="")
    return 0


def command_export_public(args: argparse.Namespace) -> int:
    repository = _repository_from_args(args)
    bundle = build_public_bundle(repository)
    output_path = Path(args.out).expanduser()
    write_public_bundle(output_path, bundle, args.format)
    print(str(output_path))
    return 0


def command_stats(args: argparse.Namespace) -> int:
    repository = _repository_from_args(args)
    metrics = collect_code_metrics(repository.repo_root)
    if args.write_metadata:
        repository.write_metrics(metrics)
    print(json.dumps(metrics.as_dict(), indent=2, sort_keys=False))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=PREFERRED_IP_CLI_NAME)
    parser.add_argument(
        "--repo-root",
        default=str(discover_repo_root()),
        help="Path to the Feel++ repository root",
    )
    parser.add_argument(
        "--metadata-path",
        default=None,
        help="Path to the public metadata YAML file. Defaults to metadata/software.public.yml.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_parser = subparsers.add_parser(
        "validate-public",
        help="Validate the public IP/APP metadata YAML file",
    )
    validate_parser.set_defaults(func=command_validate_public)

    show_parser = subparsers.add_parser(
        "show-public",
        help="Show the canonical public IP/APP metadata",
    )
    show_parser.add_argument(
        "--format",
        choices=("yaml", "json"),
        default="yaml",
        help="Output format",
    )
    show_parser.set_defaults(func=command_show_public)

    export_parser = subparsers.add_parser(
        "export-public",
        help="Export the stable public APP bundle",
    )
    export_parser.add_argument(
        "--format",
        choices=("yaml", "json"),
        required=True,
        help="Export format",
    )
    export_parser.add_argument("--out", required=True, help="Output path")
    export_parser.set_defaults(func=command_export_public)

    stats_parser = subparsers.add_parser(
        "stats",
        help="Compute best-effort public source metrics",
    )
    stats_parser.add_argument(
        "--write-metadata",
        action="store_true",
        help="Update the metrics block in metadata/software.public.yml",
    )
    stats_parser.set_defaults(func=command_stats)

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
