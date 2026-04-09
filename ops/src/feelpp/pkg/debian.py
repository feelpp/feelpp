from __future__ import annotations

from pathlib import Path
import re

from .config import PackagingContext
from .graph import build_plan, load_manifest


BUILD_DEP_FIELDS = ("Build-Depends", "Build-Depends-Indep")
PACKAGE_FIELD = "Package"
_MULTIARCH_SUFFIX_RE = re.compile(r":[A-Za-z0-9-]+$")


def control_paths_for_context(context: PackagingContext) -> list[Path]:
    manifest = load_manifest(context.manifest_path)
    plan = build_plan(manifest, dist=context.dist)
    return control_paths_for_components(context, [component.name for component in plan.components])


def control_paths_for_components(
    context: PackagingContext, components: list[str] | tuple[str, ...]
) -> list[Path]:
    paths: list[Path] = []
    for component_name in components:
        control_path = (
            context.repo_root
            / "packaging"
            / "debian"
            / component_name
            / context.dist
            / "debian"
            / "control"
        )
        if not control_path.is_file():
            raise FileNotFoundError(f"Debian control file not found: {control_path}")
        paths.append(control_path)
    return paths


def parse_control_paragraphs(path: Path) -> list[dict[str, str]]:
    paragraphs: list[dict[str, str]] = []
    current: dict[str, str] = {}
    current_field: str | None = None

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped:
            if current:
                paragraphs.append(current)
                current = {}
                current_field = None
            continue

        if raw_line.lstrip().startswith("#"):
            continue

        if raw_line[0].isspace():
            if current_field is None:
                continue
            current[current_field] = f"{current[current_field]} {stripped}"
            continue

        if ":" not in raw_line:
            continue

        key, value = raw_line.split(":", 1)
        current_field = key.strip()
        current[current_field] = value.strip()

    if current:
        paragraphs.append(current)

    return paragraphs


def _split_top_level(raw: str, separator: str) -> list[str]:
    parts: list[str] = []
    start = 0
    paren_depth = 0
    bracket_depth = 0
    angle_depth = 0

    for index, char in enumerate(raw):
        if char == "(":
            paren_depth += 1
        elif char == ")" and paren_depth:
            paren_depth -= 1
        elif char == "[":
            bracket_depth += 1
        elif char == "]" and bracket_depth:
            bracket_depth -= 1
        elif char == "<":
            angle_depth += 1
        elif char == ">" and angle_depth:
            angle_depth -= 1
        elif (
            char == separator
            and paren_depth == 0
            and bracket_depth == 0
            and angle_depth == 0
        ):
            part = raw[start:index].strip()
            if part:
                parts.append(part)
            start = index + 1

    tail = raw[start:].strip()
    if tail:
        parts.append(tail)
    return parts


def _normalize_dependency_name(raw: str) -> str | None:
    value = re.sub(r"<[^>]*>", "", raw)
    value = re.sub(r"\[[^\]]*\]", "", value)
    value = re.sub(r"\([^)]*\)", "", value)
    value = value.strip()
    if not value:
        return None
    name = value.split()[0]
    if name.startswith("${"):
        return None
    return _MULTIARCH_SUFFIX_RE.sub("", name)


def _relation_candidates(raw: str) -> list[str]:
    candidates: list[str] = []
    for alternative in _split_top_level(raw, "|"):
        normalized = _normalize_dependency_name(alternative)
        if normalized:
            candidates.append(normalized)
    return candidates


def _ordered_unique(values: list[str]) -> list[str]:
    ordered: list[str] = []
    seen: set[str] = set()
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def collect_binary_packages(control_paths: list[Path]) -> list[str]:
    binary_packages: list[str] = []
    for control_path in control_paths:
        for paragraph in parse_control_paragraphs(control_path)[1:]:
            package_name = paragraph.get(PACKAGE_FIELD)
            if package_name:
                binary_packages.append(package_name)
    return _ordered_unique(binary_packages)


def collect_seed_build_dependencies(context: PackagingContext) -> list[str]:
    control_paths = control_paths_for_context(context)
    return collect_build_dependencies_for_controls(control_paths)


def component_scope(context: PackagingContext, component: str) -> list[str]:
    manifest = load_manifest(context.manifest_path)
    ordered = list(manifest.default_components)
    if component not in ordered:
        raise ValueError(f"Unknown component: {component}")
    return ordered[: ordered.index(component) + 1]


def collect_component_build_dependencies(
    context: PackagingContext,
    component: str,
    *,
    include_prefix_scope: bool = True,
) -> list[str]:
    components = component_scope(context, component) if include_prefix_scope else [component]
    control_paths = control_paths_for_components(context, components)
    return collect_build_dependencies_for_controls(control_paths)


def collect_component_internal_build_dependencies(
    context: PackagingContext,
    component: str,
    *,
    include_prefix_scope: bool = True,
) -> list[str]:
    components = component_scope(context, component) if include_prefix_scope else [component]
    control_paths = control_paths_for_components(context, components)
    return collect_internal_build_dependencies_for_controls(control_paths)


def collect_build_dependencies_for_controls(control_paths: list[Path]) -> list[str]:
    internal_packages = set(collect_binary_packages(control_paths))
    resolved: list[str] = []

    for control_path in control_paths:
        paragraphs = parse_control_paragraphs(control_path)
        if not paragraphs:
            continue
        source_paragraph = paragraphs[0]
        for field_name in BUILD_DEP_FIELDS:
            raw_value = source_paragraph.get(field_name, "")
            if not raw_value:
                continue
            for relation in _split_top_level(raw_value, ","):
                for candidate in _relation_candidates(relation):
                    if candidate in internal_packages:
                        continue
                    resolved.append(candidate)
                    break

    return _ordered_unique(resolved)


def collect_internal_build_dependencies_for_controls(control_paths: list[Path]) -> list[str]:
    internal_packages = set(collect_binary_packages(control_paths))
    resolved: list[str] = []

    for control_path in control_paths:
        paragraphs = parse_control_paragraphs(control_path)
        if not paragraphs:
            continue
        source_paragraph = paragraphs[0]
        for field_name in BUILD_DEP_FIELDS:
            raw_value = source_paragraph.get(field_name, "")
            if not raw_value:
                continue
            for relation in _split_top_level(raw_value, ","):
                for candidate in _relation_candidates(relation):
                    if candidate not in internal_packages:
                        continue
                    resolved.append(candidate)
                    break

    return _ordered_unique(resolved)
