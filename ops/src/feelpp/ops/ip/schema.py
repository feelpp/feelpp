from __future__ import annotations

from collections.abc import Iterable
from typing import Any


SCHEMA_VERSION = 1

REQUIRED_TOP_LEVEL_KEYS = (
    "schema_version",
    "software",
    "repository",
    "languages",
    "build_tools",
    "public_license",
    "architecture",
    "public_funding",
    "public_outputs",
    "metrics",
    "private_boundary",
)

FORBIDDEN_PRIVATE_KEYS = {
    "hr",
    "inventor",
    "inventors",
    "ownership",
    "ownership_share",
    "ownership_shares",
    "payroll",
    "personal_address",
    "personal_addresses",
    "salaries",
    "salary",
    "signature",
    "signatures",
}


def _type_name(value: object) -> str:
    return type(value).__name__


def _require_mapping(payload: dict[str, Any], key: str) -> dict[str, Any]:
    value = payload.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"{key} must be a mapping, got {_type_name(value)}")
    return value


def _require_list(payload: dict[str, Any], key: str) -> list[Any]:
    value = payload.get(key)
    if not isinstance(value, list):
        raise ValueError(f"{key} must be a list, got {_type_name(value)}")
    return value


def _require_string(payload: dict[str, Any], key: str, *, allow_null: bool = False) -> None:
    value = payload.get(key)
    if allow_null and value is None:
        return
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{key} must be a non-empty string")


def _require_string_list(payload: dict[str, Any], key: str) -> None:
    values = _require_list(payload, key)
    for index, value in enumerate(values):
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{key}[{index}] must be a non-empty string")


def _iter_mapping_keys(value: Any) -> Iterable[str]:
    if isinstance(value, dict):
        for key, child in value.items():
            yield str(key)
            yield from _iter_mapping_keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from _iter_mapping_keys(child)


def _normalize_key(key: str) -> str:
    return key.strip().lower().replace("-", "_")


def validate_public_metadata(payload: dict[str, Any]) -> None:
    missing = [key for key in REQUIRED_TOP_LEVEL_KEYS if key not in payload]
    if missing:
        raise ValueError(f"Missing required top-level key(s): {', '.join(missing)}")

    if payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"schema_version must be {SCHEMA_VERSION}")

    for key in _iter_mapping_keys(payload):
        normalized = _normalize_key(key)
        if normalized in FORBIDDEN_PRIVATE_KEYS:
            raise ValueError(f"Private metadata key is not allowed in public metadata: {key}")

    software = _require_mapping(payload, "software")
    _require_string(software, "name")
    _require_string(software, "description")
    _require_string(software, "homepage", allow_null=True)

    repository = _require_mapping(payload, "repository")
    for key in ("hosting", "owner", "name", "url", "default_branch"):
        _require_string(repository, key)

    languages = _require_list(payload, "languages")
    for index, language in enumerate(languages):
        if not isinstance(language, dict):
            raise ValueError(f"languages[{index}] must be a mapping")
        _require_string(language, "name")

    build_tools = _require_list(payload, "build_tools")
    for index, tool in enumerate(build_tools):
        if not isinstance(tool, dict):
            raise ValueError(f"build_tools[{index}] must be a mapping")
        _require_string(tool, "name")

    public_license = _require_mapping(payload, "public_license")
    _require_string_list(public_license, "expressions")
    files = _require_list(public_license, "files")
    for index, file_entry in enumerate(files):
        if not isinstance(file_entry, dict):
            raise ValueError(f"public_license.files[{index}] must be a mapping")
        _require_string(file_entry, "path")

    architecture = _require_mapping(payload, "architecture")
    _require_string(architecture, "overview")
    components = _require_list(architecture, "components")
    for index, component in enumerate(components):
        if not isinstance(component, dict):
            raise ValueError(f"architecture.components[{index}] must be a mapping")
        _require_string(component, "name")
        _require_string(component, "description")

    public_funding = _require_mapping(payload, "public_funding")
    _require_list(public_funding, "grants")

    public_outputs = _require_mapping(payload, "public_outputs")
    for key, value in public_outputs.items():
        if value is not None and not isinstance(value, dict):
            raise ValueError(f"public_outputs.{key} must be a mapping or null")

    metrics = _require_mapping(payload, "metrics")
    for key in ("approximate_lines_of_code", "approximate_bytes", "counted_files"):
        value = metrics.get(key)
        if value is not None and (not isinstance(value, int) or value < 0):
            raise ValueError(f"metrics.{key} must be a non-negative integer or null")

    private_boundary = _require_mapping(payload, "private_boundary")
    _require_string(private_boundary, "private_dossier_repository")
    _require_string_list(private_boundary, "excluded_information")
