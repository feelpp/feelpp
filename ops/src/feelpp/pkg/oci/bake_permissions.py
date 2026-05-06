from __future__ import annotations

import json
from pathlib import Path


def _is_non_filesystem_reference(value: str) -> bool:
    raw = value.strip()
    return (
        not raw
        or raw.startswith("target:")
        or "://" in raw
    )


def _iter_target_path_values(target_spec: dict[str, object]) -> list[str]:
    values: list[str] = []
    for key in ("context", "dockerfile"):
        value = target_spec.get(key)
        if isinstance(value, str):
            values.append(value)

    contexts = target_spec.get("contexts")
    if isinstance(contexts, dict):
        for value in contexts.values():
            if isinstance(value, str):
                values.append(value)

    return values


def _selected_target_specs(
    payload: dict[str, object],
    *,
    groups: list[str] | None,
) -> list[dict[str, object]]:
    targets = payload.get("target", {})
    if not isinstance(targets, dict):
        return []

    if not groups:
        return [spec for spec in targets.values() if isinstance(spec, dict)]

    raw_groups = payload.get("group", {})
    group_map = raw_groups if isinstance(raw_groups, dict) else {}
    selected_specs: list[dict[str, object]] = []
    seen_targets: set[str] = set()
    seen_groups: set[str] = set()

    def visit(name: str) -> None:
        if name in targets:
            if name in seen_targets:
                return
            seen_targets.add(name)
            target_spec = targets[name]
            if isinstance(target_spec, dict):
                selected_specs.append(target_spec)
            return

        if name in seen_groups:
            return
        seen_groups.add(name)
        group_spec = group_map.get(name)
        if not isinstance(group_spec, dict):
            return
        group_targets = group_spec.get("targets", [])
        if not isinstance(group_targets, list):
            return
        for target_name in group_targets:
            if isinstance(target_name, str):
                visit(target_name)

    for group_name in groups:
        visit(group_name)

    return selected_specs


def required_fs_read_paths_from_bake_payload(
    payload: dict[str, object],
    *,
    base_dir: Path,
    groups: list[str] | None = None,
) -> list[str]:
    resolved_base_dir = base_dir.resolve()
    required_paths: list[str] = []
    for target_spec in _selected_target_specs(payload, groups=groups):
        for raw_value in _iter_target_path_values(target_spec):
            if _is_non_filesystem_reference(raw_value):
                continue

            candidate = Path(raw_value)
            if candidate.is_absolute():
                required_paths.append(str(candidate.resolve()))
                continue

            resolved_candidate = (resolved_base_dir / candidate).resolve()
            try:
                resolved_candidate.relative_to(resolved_base_dir)
            except ValueError:
                required_paths.append(str(resolved_candidate))

    return sorted(set(required_paths))


def required_fs_read_paths_from_bake_file(
    bake_file: Path,
    *,
    groups: list[str] | None = None,
) -> list[str]:
    payload = json.loads(bake_file.read_text(encoding="utf-8"))
    return required_fs_read_paths_from_bake_payload(payload, base_dir=bake_file.parent, groups=groups)


def fs_read_allow_flags(paths: list[str]) -> list[str]:
    return [f"--allow=fs.read={path}" for path in sorted(set(paths))]
