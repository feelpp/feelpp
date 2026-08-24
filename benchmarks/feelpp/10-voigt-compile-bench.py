#!/usr/bin/env python3

import argparse
import json
import pathlib
import shlex
import subprocess
import sys
import time


def load_compile_template(build_dir: pathlib.Path) -> tuple[pathlib.Path, list[str]]:
    compile_commands = build_dir / "compile_commands.json"
    if not compile_commands.exists():
        raise FileNotFoundError(f"missing {compile_commands}")

    entries = json.loads(compile_commands.read_text())
    for entry in entries:
        source = pathlib.Path(entry["file"])
        if source.as_posix().endswith("benchmarks/feelpp/08-voigt-elasticity.cpp"):
            if "arguments" in entry:
                return pathlib.Path(entry["directory"]), list(entry["arguments"])
            return pathlib.Path(entry["directory"]), shlex.split(entry["command"])
    raise RuntimeError("could not find a compile_commands entry for benchmarks/feelpp/08-voigt-elasticity.cpp")


def rewritten_compile_command(template: list[str], source: pathlib.Path, output: pathlib.Path) -> list[str]:
    command = list(template)
    has_compile_only = False
    for idx, token in enumerate(command):
        if token == "-c" and idx + 1 < len(command):
            command[idx + 1] = str(source)
            has_compile_only = True
        elif token == "-o" and idx + 1 < len(command):
            command[idx + 1] = str(output)
    if not has_compile_only:
        command.extend(["-c", str(source)])
    command.extend(["-fsyntax-only"])
    return command


def main() -> int:
    parser = argparse.ArgumentParser(description="Measure compile time for representative Voigt DSEL translation units.")
    parser.add_argument("build_dir", nargs="?", default="build/default")
    args = parser.parse_args()

    repo_root = pathlib.Path(__file__).resolve().parents[2]
    build_dir = (repo_root / args.build_dir).resolve()
    compile_dir = pathlib.Path(__file__).resolve().parent / "compile"
    out_dir = build_dir / "benchmarks" / "feelpp" / "compile_bench"
    out_dir.mkdir(parents=True, exist_ok=True)

    workdir, template = load_compile_template(build_dir)

    sources = [
        compile_dir / "tensor_elasticity_expr.cpp",
        compile_dir / "mandel_elasticity_expr.cpp",
        compile_dir / "voigt_elasticity_expr.cpp",
        compile_dir / "inverse_conversion_expr.cpp",
    ]

    print("Voigt compile benchmark", flush=True)
    print(f"build_dir: {build_dir}", flush=True)
    for source in sources:
        output = out_dir / f"{source.stem}.o"
        if output.exists():
            output.unlink()

        command = rewritten_compile_command(template, source, output)
        started = time.perf_counter()
        result = subprocess.run(command, cwd=workdir, capture_output=True, text=True)
        elapsed = time.perf_counter() - started

        if result.returncode != 0:
            sys.stderr.write(result.stdout)
            sys.stderr.write(result.stderr)
            raise SystemExit(result.returncode)

        print(f"{source.stem}: {elapsed:.3f}s", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
