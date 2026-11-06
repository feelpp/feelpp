#!/usr/bin/env python3
# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Solve the curved P2/G2 affine HDG case and check the actual EnSight fields.

VTK is optional: return 77 when its Python reader is unavailable. Temporary
solver output follows Python's tempfile policy, including the TMPDIR setting.
"""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile


def check_output(case_file, reader_type, tolerance):
    reader = reader_type()
    errors = []
    reader.AddObserver("ErrorEvent", lambda *args: errors.append("VTK reader error"))
    reader.SetCaseFileName(str(case_file))
    reader.ReadAllVariablesOn()
    reader.Update()
    if errors:
        raise RuntimeError(f"could not read {case_file}: {errors}")

    iterator = reader.GetOutput().NewIterator()
    iterator.InitTraversal()
    total_cells = 0
    max_error = {name: 0.0 for name in ("potential", "potential.exact", "potentialpp", "flux")}
    max_midpoint_displacement = 0.0
    while not iterator.IsDoneWithTraversal():
        block = iterator.GetCurrentDataObject()
        iterator.GoToNextItem()
        if block is None or not hasattr(block, "GetNumberOfCells"):
            continue
        cells = block.GetNumberOfCells()
        if cells == 0:
            continue
        data = block.GetCellData()
        fields = {name: data.GetArray(name) for name in max_error}
        for name, field in fields.items():
            components = 3 if name == "flux" else 1
            if field is None or field.GetNumberOfTuples() != cells or field.GetNumberOfComponents() != components:
                raise RuntimeError(f"missing or invalid cell field {name!r}")

        for index in range(cells):
            cell = block.GetCell(index)
            # VTK_QUADRATIC_TRIANGLE: vertices first, then edge midpoints.
            if cell.GetCellType() != 22 or cell.GetNumberOfPoints() != 6:
                raise RuntimeError("expected quadratic triangles from geometry.order=2")
            ids = [cell.GetPointId(j) for j in range(6)]
            if len(set(ids)) != 6 or any(not 0 <= node < block.GetNumberOfPoints() for node in ids):
                raise RuntimeError("invalid quadratic connectivity")
            points = [block.GetPoint(node) for node in ids]
            if not all(math.isfinite(x) for point in points for x in point):
                raise RuntimeError("nonfinite exported geometry")
            for a, b, middle in ((0, 1, 3), (1, 2, 4), (2, 0, 5)):
                length = math.dist(points[a], points[b])
                if length == 0:
                    raise RuntimeError("degenerate edge in exported geometry")
                midpoint = [(x + y) / 2 for x, y in zip(points[a], points[b])]
                displacement = math.dist(points[middle], midpoint) / length
                max_midpoint_displacement = max(max_midpoint_displacement, displacement)

            # P0 export samples the reference centroid. For a quadratic map,
            # its physical image is not the arithmetic mean of the vertices.
            weights = (-1 / 9, -1 / 9, -1 / 9, 4 / 9, 4 / 9, 4 / 9)
            expected = sum(weight * (point[0] + point[1]) for weight, point in zip(weights, points))
            for name in ("potential", "potential.exact", "potentialpp"):
                value = fields[name].GetTuple1(index)
                if not math.isfinite(value):
                    raise RuntimeError(f"nonfinite {name} at cell {index}")
                max_error[name] = max(max_error[name], abs(value - expected))
            for value, expected_flux in zip(fields["flux"].GetTuple(index), (-1, -1, 0)):
                if not math.isfinite(value):
                    raise RuntimeError(f"nonfinite flux at cell {index}")
                max_error["flux"] = max(max_error["flux"], abs(value - expected_flux))
            if abs(fields["potentialpp"].GetTuple1(index) - fields["potential.exact"].GetTuple1(index)) > tolerance:
                raise RuntimeError(f"exported potentialpp differs from potential.exact at cell {index}")
        total_cells += cells

    if total_cells == 0:
        raise RuntimeError("no cells were read from the exported case")
    if not 1e-4 < max_midpoint_displacement < 0.2:
        raise RuntimeError(f"expected a curved circle mesh with consistent connectivity; midpoint ratio={max_midpoint_displacement}")
    if any(error > tolerance for error in max_error.values()):
        raise RuntimeError(f"exported affine-field errors exceed {tolerance}: {max_error}")
    return {"cells": total_cells, "max_error": max_error, "max_midpoint_ratio": max_midpoint_displacement}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--checker", type=Path, required=True, help="checker JSON file")
    parser.add_argument("--geometry", type=Path, required=True, help="curved circle .geo file")
    parser.add_argument("--checker-script", type=Path, default=Path(__file__).resolve().parent / "python" / "laplacian.py")
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--tolerance", type=float, default=2e-6, help="absolute tolerance for EnSight float values")
    args = parser.parse_args()
    if not math.isfinite(args.tolerance) or args.tolerance <= 0 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("tolerance and timeout must be positive finite values")
    try:
        from vtkmodules.vtkIOEnSight import vtkEnSightGoldBinaryReader
    except ImportError:
        print("SKIP: VTK Python EnSight reader is unavailable")
        return 77

    paths = {name: getattr(args, name).resolve() for name in ("executable", "config", "checker", "geometry", "checker_script")}
    for name, path in paths.items():
        if not path.is_file():
            parser.error(f"{name} does not exist: {path}")
    with tempfile.TemporaryDirectory(prefix="feelpp-hdg-export-") as temporary:
        directory = Path(temporary)
        command = [
            str(paths["executable"]), f"--config-file={paths['config']}",
            f"--directory={directory / 'results'}", f"--gmsh.filename={paths['geometry']}",
            "--order=2", "--geometry.order=2", "--gmsh.hsize=0.4", "--sc.condense=true",
            "--k=1", "--solution.sympy.p=x+y", "--solution.p=x+y:x:y", "--solution.u={-1,-1}",
            "--functions.f=0", "--checker.solution=x+y", "--checker.check=true",
            "--checker.tolerance.exact=1e-9", f"--checker.filename={paths['checker']}",
            f"--checker.script={paths['checker_script']}", "--exporter.format=ensightgold",
            "--exporter.geometry=static", "--exporter.export=true",
        ]
        completed = subprocess.run(command, cwd=directory, text=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, timeout=args.timeout, check=False)
        try:
            if completed.returncode:
                raise RuntimeError(f"HDG solver exited with {completed.returncode}")
            cases = list(directory.rglob("*.case"))
            if len(cases) != 1:
                raise RuntimeError(f"expected one fresh EnSight case, found {cases}")
            result = check_output(cases[0], vtkEnSightGoldBinaryReader, args.tolerance)
        except Exception:
            sys.stderr.write(completed.stdout)
            raise
        print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"HDG export check failed: {error}", file=sys.stderr)
        sys.exit(1)
