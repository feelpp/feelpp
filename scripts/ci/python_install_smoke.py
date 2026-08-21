#!/usr/bin/env python3

import argparse
import importlib
import pathlib
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description="Import installed Feel++ Python modules from a prefix")
    parser.add_argument("--prefix", required=True, help="Installation prefix used for the smoke run")
    parser.add_argument("--module", action="append", required=True, dest="modules", help="Module to import")
    args = parser.parse_args()

    prefix = pathlib.Path(args.prefix).resolve()

    for module_name in args.modules:
        module = importlib.import_module(module_name)
        module_file = getattr(module, "__file__", None)
        if not module_file:
            print(f"{module_name}: imported")
            continue

        module_path = pathlib.Path(module_file).resolve()
        if not module_path.is_relative_to(prefix):
            print(f"{module_name}: imported from unexpected location {module_path}", file=sys.stderr)
            return 1

        print(f"{module_name}: {module_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
