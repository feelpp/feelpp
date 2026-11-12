#!/usr/bin/env python3

# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#
# SPDX-FileCopyrightText: 2026 University of Strasbourg
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Reject runtime-named public callables that duplicate semantic APIs."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


HEADER_SUFFIXES = {".h", ".hh", ".hpp", ".hxx"}
CALLABLE = re.compile(r"\b([A-Za-z_][A-Za-z0-9_]*)\s*(?=\()")

# These names describe distinct runtime-sized algorithms, runtime capability
# queries, or the validated construction tag. They are not alternate semantic
# accessors for order, cardinality, or DOF topology.
ALLOWED_RUNTIME_CALLABLES: dict[str, frozenset[str] | None] = {
    "RuntimeOrder": None,
    "basisDerivateRuntimeLowOrderDispatch": frozenset(
        {"feelpp/feel/feelpoly/orthonormalpolynomialset.hpp"}
    ),
    "basisEvaluateRuntimeLowOrderDispatch": frozenset(
        {"feelpp/feel/feelpoly/orthonormalpolynomialset.hpp"}
    ),
    "derivateRuntime": frozenset(
        {"feelpp/feel/feelpoly/dubiner.hpp", "feelpp/feel/feelpoly/legendre.hpp"}
    ),
    "evaluateRuntime": frozenset(
        {"feelpp/feel/feelpoly/dubiner.hpp", "feelpp/feel/feelpoly/legendre.hpp"}
    ),
    "petscHasHpddmAuxiliaryMatRuntime": frozenset(
        {"feelpp/feel/feelalg/petschpddm.hpp"}
    ),
    "petscHasHpddmRuntime": frozenset(
        {"feelpp/feel/feelalg/petschpddm.hpp"}
    ),
    "scalingsRuntime": frozenset(
        {
            "feelpp/feel/feelpoly/dubiner.hpp",
            "feelpp/feel/feelpoly/expansions.hpp",
        }
    ),
}


def strip_comments_and_literals(source: str) -> str:
    """Replace comments and literals with spaces while preserving newlines."""
    result = list(source)
    index = 0
    size = len(source)

    while index < size:
        if source.startswith("//", index):
            end = source.find("\n", index + 2)
            if end == -1:
                end = size
            result[index:end] = " " * (end - index)
            index = end
        elif source.startswith("/*", index):
            end = source.find("*/", index + 2)
            end = size if end == -1 else end + 2
            for offset in range(index, end):
                if source[offset] != "\n":
                    result[offset] = " "
            index = end
        elif source[index] in {'"', "'"}:
            quote = source[index]
            result[index] = " "
            index += 1
            while index < size:
                if source[index] == "\\":
                    result[index] = " "
                    if index + 1 < size:
                        if source[index + 1] != "\n":
                            result[index + 1] = " "
                        index += 2
                    else:
                        index += 1
                elif source[index] == quote:
                    result[index] = " "
                    index += 1
                    break
                else:
                    if source[index] != "\n":
                        result[index] = " "
                    index += 1
        else:
            index += 1

    return "".join(result)


def forbidden_callables(header: Path, source_root: Path) -> list[str]:
    source = strip_comments_and_literals(
        header.read_text(encoding="utf-8", errors="replace")
    )
    relative = header.relative_to(source_root).as_posix()
    diagnostics: list[str] = []

    for match in CALLABLE.finditer(source):
        name = match.group(1)
        if "runtime" not in name.lower():
            continue
        if name.startswith("M_"):
            continue
        if source[max(0, match.start() - 5) : match.start()] == "std::":
            continue
        allowed_paths = ALLOWED_RUNTIME_CALLABLES.get(name, frozenset())
        if allowed_paths is None or relative in allowed_paths:
            continue

        line = source.count("\n", 0, match.start()) + 1
        diagnostics.append(
            f"{relative}:{line}: forbidden runtime-named callable '{name}'"
        )

    return diagnostics


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "source_root", type=Path, help="Feel++ repository source directory"
    )
    args = parser.parse_args()

    source_root = args.source_root.resolve()
    public_root = source_root / "feelpp" / "feel"
    headers = sorted(
        path
        for path in public_root.rglob("*")
        if path.is_file() and path.suffix in HEADER_SUFFIXES
    )

    diagnostics = [
        diagnostic
        for header in headers
        for diagnostic in forbidden_callables(header, source_root)
    ]
    if diagnostics:
        print("Runtime API naming check failed:")
        print("\n".join(diagnostics))
        print(
            "Use the common semantic accessor name, or document a genuinely "
            "distinct runtime algorithm in ALLOWED_RUNTIME_CALLABLES."
        )
        return 1

    print(f"Runtime API naming check passed for {len(headers)} public headers.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
