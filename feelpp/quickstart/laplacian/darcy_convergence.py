#!/usr/bin/env python3

# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#
# SPDX-FileCopyrightText: 2026 University of Strasbourg
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import argparse
import csv
import hashlib
import json
import math
import pathlib
import re
import subprocess
import sys


METRICS = {
    "flux_l2": r"\|\s*\|\|j-j_h\|\|_L2\s*\|\s*([0-9.eE+-]+)",
    "div_l2": r"\|\s*\|\|div\(j_h\)-g\|\|_L2\s*\|\s*([0-9.eE+-]+)",
    "pressure_l2": r"\|\s*\|\|V-V_h\|\|_L2\s*\|\s*([0-9.eE+-]+)",
    "mass_residual": r"\|\s*weak mass residual \(Euclidean norm\)\s*\|\s*([0-9.eE+-]+)",
    "boundary_flux_residual": r"\|\s*\|\|j_h\.n-j_N\|\|_L2\(Neumann\)\s*\|\s*([0-9.eE+-]+)",
}


def parse_pair(value):
    family, order = value.split( ":", 1 )
    if family not in ("rt", "bdm"):
        raise argparse.ArgumentTypeError(f"unsupported family {family}")
    return family, int(order)


def parse_metrics(output):
    result = {}
    for name, pattern in METRICS.items():
        match = re.search(pattern, output)
        if not match:
            raise RuntimeError(f"metric {name} was not found in solver output")
        result[name] = float(match.group(1))
    return result


def run_case(executable, config, run_id, formulation, family, order, hsize):
    directory = f"qs_darcy_convergence/{run_id}/{formulation}/{family}{order}/h{hsize:g}"
    command = [
        str(executable),
        "--config-file", str(config),
        f"--directory={directory}",
        f"--gmsh.hsize={hsize}",
        f"--darcy.flux-family={family}",
        f"--darcy.order={order}",
        "--darcy.convergence-case=true",
        "--ksp-monitor=false",
        "--no-export",
    ]
    if formulation.startswith("hybrid"):
        command.append(f"--sc.condense={'true' if formulation == 'hybrid-condensed' else 'false'}")
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise RuntimeError(
            f"{formulation} {family.upper()}{order} h={hsize:g} failed with exit code {completed.returncode}"
        )
    return parse_metrics(completed.stdout)


def expected_rates(family, order):
    if family == "rt":
        return {"flux_l2": order + 1, "div_l2": order + 1, "pressure_l2": order + 1}
    return {"flux_l2": order + 1, "div_l2": order, "pressure_l2": order}


def add_rates(rows):
    grouped = {}
    for row in rows:
        grouped.setdefault((row["formulation"], row["family"], row["order"]), []).append(row)
    for values in grouped.values():
        values.sort(key=lambda row: row["h"], reverse=True)
        previous = None
        for row in values:
            for metric in ("flux_l2", "div_l2", "pressure_l2"):
                row[f"{metric}_rate"] = ""
            if previous is not None:
                denominator = math.log(previous["h"] / row["h"])
                for metric in ("flux_l2", "div_l2", "pressure_l2"):
                    if row[metric] > 0 and previous[metric] > 0:
                        row[f"{metric}_rate"] = math.log(previous[metric] / row[metric]) / denominator
            previous = row


def check_results(rows, rate_slack, comparison_tolerance):
    failures = []
    grouped = {}
    for row in rows:
        grouped.setdefault((row["formulation"], row["family"], row["order"]), []).append(row)
        if not math.isfinite(row["mass_residual"]) or row["mass_residual"] > 1e-8:
            failures.append(
                f"{row['formulation']} {row['family'].upper()}{row['order']} h={row['h']}: "
                f"weak mass residual must be finite and <= 1e-8; got {row['mass_residual']:.6e}"
            )

    for key, values in grouped.items():
        values.sort(key=lambda row: row["h"], reverse=True)
        expected = expected_rates(key[1], key[2])
        if len(values) < 3:
            failures.append(f"{key}: at least three refinement levels are required")
            continue
        for metric, expected_rate in expected.items():
            rates = [row[f"{metric}_rate"] for row in values if row[f"{metric}_rate"] != ""]
            usable = [rate for rate in rates[-2:] if math.isfinite(rate)]
            if not usable or min(usable) < expected_rate - rate_slack:
                failures.append(
                    f"{key[0]} {key[1].upper()}{key[2]} {metric}: expected >= "
                    f"{expected_rate-rate_slack:.3f}, observed {usable}"
                )

    by_discretization = {}
    for row in rows:
        by_discretization.setdefault((row["family"], row["order"], row["h"]), []).append(row)
    for key, values in by_discretization.items():
        reference = next((row for row in values if row["formulation"] == "conforming"), None)
        if reference is None:
            continue
        for row in values:
            if row is reference:
                continue
            for metric in ("flux_l2", "div_l2", "pressure_l2"):
                scale = max(reference[metric], row[metric], 1e-14)
                relative = abs(reference[metric] - row[metric]) / scale
                if relative > comparison_tolerance:
                    failures.append(
                        f"{key} {metric}: conforming={reference[metric]:.6e}, "
                        f"{row['formulation']}={row[metric]:.6e}, relative difference={relative:.3e}"
                    )
    return failures


def write_plotly_html(rows, output):
    metrics = (
        ("flux_l2", "Flux L2 error"),
        ("div_l2", "Divergence L2 error"),
        ("pressure_l2", "Pressure L2 error"),
    )
    grouped = {}
    for row in rows:
        grouped.setdefault((row["formulation"], row["family"], row["order"]), []).append(row)

    plots = []
    for metric, title in metrics:
        traces = []
        for (formulation, family, order), values in sorted(grouped.items()):
            values.sort(key=lambda row: row["h"])
            expected = expected_rates(family, order)[metric]
            traces.append({
                "x": [row["h"] for row in values],
                "y": [row[metric] for row in values],
                "mode": "lines+markers",
                "name": f"{formulation} {family.upper()}{order} (expected {expected})",
                "hovertemplate": "h=%{x:.4g}<br>error=%{y:.6e}<extra>%{fullData.name}</extra>",
            })
        plots.append({"metric": metric, "title": title, "traces": traces})

    output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(plots)
    output.write_text(f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <title>Darcy RT/BDM convergence</title>
  <script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>
  <style>
    body {{ font-family: sans-serif; margin: 1.5rem; }}
    .plot {{ width: 100%; height: 34rem; margin-bottom: 2rem; }}
  </style>
</head>
<body>
  <h1>Darcy RT/BDM convergence</h1>
  <p>Interactive log-log plots generated from the Phase 2b CSV data.</p>
  <div id=\"flux_l2\" class=\"plot\"></div>
  <div id=\"div_l2\" class=\"plot\"></div>
  <div id=\"pressure_l2\" class=\"plot\"></div>
  <script>
    const plots = {payload};
    for (const plot of plots) {{
      Plotly.newPlot(plot.metric, plot.traces, {{
        title: plot.title,
        xaxis: {{ title: 'requested h', type: 'log', autorange: 'reversed' }},
        yaxis: {{ title: 'L2 error', type: 'log' }},
        legend: {{ orientation: 'h', y: -0.25 }},
        margin: {{ b: 150 }}
      }}, {{ responsive: true }});
    }}
  </script>
</body>
</html>
""", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Darcy RT/BDM convergence study")
    parser.add_argument("--conforming", type=pathlib.Path, required=True)
    parser.add_argument("--hybrid", type=pathlib.Path, required=True)
    parser.add_argument("--config", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--plot", type=pathlib.Path)
    parser.add_argument("--pair", action="append", type=parse_pair, required=True)
    parser.add_argument("--h", action="append", type=float, required=True)
    parser.add_argument(
        "--formulation", action="append",
        choices=("conforming", "hybrid-monolithic", "hybrid-condensed"), required=True
    )
    parser.add_argument("--rate-slack", type=float, default=0.25)
    parser.add_argument("--comparison-tolerance", type=float, default=5e-5)
    parser.add_argument("--no-check", action="store_true")
    args = parser.parse_args()

    rows = []
    run_id = hashlib.sha256(args.config.read_bytes() + b"darcy-poly6-v1").hexdigest()[:12]
    for family, order in args.pair:
        for formulation in args.formulation:
            if family == "bdm" and order == 3 and formulation == "hybrid-monolithic":
                continue
            executable = args.conforming if formulation == "conforming" else args.hybrid
            for hsize in sorted(args.h, reverse=True):
                metrics = run_case(executable, args.config, run_id, formulation, family, order, hsize)
                rows.append({
                    "formulation": formulation,
                    "family": family,
                    "order": order,
                    "h": hsize,
                    **metrics,
                })
                print(
                    f"{formulation:18s} {family.upper()}{order} h={hsize:g} "
                    f"flux={metrics['flux_l2']:.6e} div={metrics['div_l2']:.6e} "
                    f"pressure={metrics['pressure_l2']:.6e}"
                )

    add_rates(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "formulation", "family", "order", "h",
        "flux_l2", "flux_l2_rate", "div_l2", "div_l2_rate",
        "pressure_l2", "pressure_l2_rate", "mass_residual", "boundary_flux_residual",
    ]
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    if args.plot:
        write_plotly_html(rows, args.plot)

    if not args.no_check:
        failures = check_results(rows, args.rate_slack, args.comparison_tolerance)
        if failures:
            for failure in failures:
                print(f"FAIL: {failure}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
