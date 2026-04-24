from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .planner_directives import write_github_output


COMPONENT_JOBS = ("feelpp", "testsuite", "toolboxes", "mor")


def _load_json_value(raw: str, *, default: object) -> object:
    if not raw.strip():
        return default
    return json.loads(raw)


def _lower_unique(values: list[str]) -> list[str]:
    seen: set[str] = set()
    normalized: list[str] = []
    for value in values:
        lowered = str(value).lower()
        if lowered in seen:
            continue
        seen.add(lowered)
        normalized.append(lowered)
    return normalized


def _as_bool(value: object, *, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"1", "true", "yes", "on"}:
            return True
        if lowered in {"0", "false", "no", "off"}:
            return False
    return default


def _resolve_ci_profile(config: dict[str, object]) -> tuple[dict[str, object], dict[str, dict[str, object]], str, str]:
    profiles = config.get("profiles")
    if not isinstance(profiles, dict):
        raise ValueError("plan-ci.json must define profiles.ci and profiles.<matrixCatalogProfile>")

    ci_profile = profiles.get("ci")
    if not isinstance(ci_profile, dict):
        raise ValueError("plan-ci.json must define profiles.ci")

    catalog_profile_name = str(ci_profile.get("matrixCatalogProfile", "")).strip().lower()
    if not catalog_profile_name:
        raise ValueError("profiles.ci.matrixCatalogProfile is required")

    catalog_profile = profiles.get(catalog_profile_name)
    if not isinstance(catalog_profile, dict):
        raise ValueError(f'profiles.{catalog_profile_name} is required')

    raw_catalog = catalog_profile.get("catalog")
    if not isinstance(raw_catalog, dict):
        raise ValueError(f'profiles.{catalog_profile_name}.catalog is required')

    catalog = {str(key).lower(): value for key, value in raw_catalog.items() if isinstance(value, dict)}
    full_build = ci_profile.get("fullBuild")
    if not isinstance(full_build, dict):
        raise ValueError("profiles.ci.fullBuild is required")

    full_job = str(full_build.get("job", "")).strip().lower()
    if not full_job:
        raise ValueError("profiles.ci.fullBuild.job is required")

    return ci_profile, catalog, catalog_profile_name, full_job


def _build_matrix(targets: list[str], catalog: dict[str, dict[str, object]], warnings: list[str]) -> dict[str, list[dict[str, object]]]:
    include: list[dict[str, object]] = []
    for target in targets:
        row = catalog.get(target)
        if row is None:
            warnings.append(f'target "{target}" is missing from the matrix catalog')
            continue
        include.append({"target": target, **row})
    return {"include": include}


def compute_workflow_plan(
    *,
    config: dict[str, object],
    mode: str,
    targets: list[str],
    enabled_jobs: list[str],
) -> dict[str, object]:
    ci_profile, catalog, _catalog_profile_name, full_job = _resolve_ci_profile(config)
    requested_mode = mode.strip().lower() or "components"
    requested_targets = _lower_unique(targets)
    requested_jobs = _lower_unique(enabled_jobs)
    requested_component_jobs = [job for job in requested_jobs if job in COMPONENT_JOBS]
    requested_full = requested_mode == "full" or full_job in requested_jobs

    raw_full_targets = ci_profile.get("fullBuild", {}).get("targets", [])
    configured_full_targets = _lower_unique(
        [target for target in raw_full_targets if isinstance(target, str)]
    )

    component_targets: list[str] = []
    full_targets: list[str] = []
    rerouted_full_targets: list[str] = []
    warnings: list[str] = []

    for target in requested_targets:
        row = catalog.get(target)
        if row is None:
            warnings.append(f'ignoring unknown target "{target}"')
            continue

        supports_components = _as_bool(
            row.get("ci_components"),
            default=str(row.get("image_backend", "")).lower() != "spack",
        )
        supports_full = _as_bool(
            row.get("ci_full"),
            default=target in configured_full_targets or str(row.get("image_backend", "")).lower() == "spack",
        )

        if requested_full:
            if supports_full:
                full_targets.append(target)
            else:
                warnings.append(f'target "{target}" does not support full builds and was skipped')
            continue

        if supports_components:
            component_targets.append(target)
            continue

        if supports_full:
            full_targets.append(target)
            rerouted_full_targets.append(target)
            continue

        warnings.append(f'target "{target}" does not support component or full CI and was skipped')

    component_targets = _lower_unique(component_targets)
    full_targets = _lower_unique(full_targets)
    rerouted_full_targets = _lower_unique(rerouted_full_targets)

    component_targets_present = bool(component_targets)
    full_targets_present = bool(full_targets)

    run_feelpp = component_targets_present and bool(
        {"feelpp", "testsuite", "toolboxes", "mor"} & set(requested_component_jobs)
    )
    run_testsuite = component_targets_present and "testsuite" in requested_component_jobs
    run_toolboxes = component_targets_present and bool({"toolboxes", "mor"} & set(requested_component_jobs))
    run_mor = component_targets_present and "mor" in requested_component_jobs
    run_full = full_targets_present and (requested_full or bool(rerouted_full_targets))

    component_matrix = _build_matrix(component_targets, catalog, warnings)
    full_matrix = _build_matrix(full_targets, catalog, warnings)

    return {
        "requested_mode": requested_mode,
        "requested_targets_json": json.dumps(requested_targets),
        "component_targets_json": json.dumps(component_targets),
        "full_targets_json": json.dumps(full_targets),
        "component_matrix_json": json.dumps(component_matrix),
        "full_matrix_json": json.dumps(full_matrix),
        "requested_component_jobs_json": json.dumps(requested_component_jobs),
        "rerouted_full_targets_json": json.dumps(rerouted_full_targets),
        "run_feelpp": str(run_feelpp).lower(),
        "run_testsuite": str(run_testsuite).lower(),
        "run_toolboxes": str(run_toolboxes).lower(),
        "run_mor": str(run_mor).lower(),
        "run_full": str(run_full).lower(),
        "component_target_count": str(len(component_matrix["include"])),
        "full_target_count": str(len(full_matrix["include"])),
        "warnings_json": json.dumps(_lower_unique(warnings)),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m feelpp.ops.ci.workflow_plan")
    parser.add_argument("--config-path", required=True, help="Path to .github/plan-ci.json")
    parser.add_argument("--mode", default="", help="Resolved planner mode")
    parser.add_argument("--targets-json", default="[]", help="Selected targets JSON array")
    parser.add_argument("--enabled-jobs-json", default="[]", help="Enabled jobs JSON array")
    parser.add_argument(
        "--github-output",
        default="",
        help="Path to the GitHub Actions output file. Defaults to $GITHUB_OUTPUT.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    config = json.loads(Path(args.config_path).read_text(encoding="utf-8"))
    targets = _load_json_value(args.targets_json, default=[])
    enabled_jobs = _load_json_value(args.enabled_jobs_json, default=[])

    if not isinstance(targets, list):
        print("--targets-json must decode to a JSON array", file=sys.stderr)
        return 1
    if not isinstance(enabled_jobs, list):
        print("--enabled-jobs-json must decode to a JSON array", file=sys.stderr)
        return 1

    try:
        outputs = compute_workflow_plan(
            config=config,
            mode=args.mode,
            targets=[str(target) for target in targets],
            enabled_jobs=[str(job) for job in enabled_jobs],
        )
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    output_path = args.github_output or ""
    if output_path:
        for name, value in outputs.items():
            write_github_output(output_path=output_path, output_name=name, value=value)
    else:
        for name, value in outputs.items():
            print(f"{name}={value}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
