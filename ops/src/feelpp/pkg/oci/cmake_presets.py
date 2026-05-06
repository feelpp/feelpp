from __future__ import annotations

from .catalog import ImageTarget


DEFAULT_COMPONENT_PRESETS = {
    "feelpp": "feelpp",
    "testsuite": "testsuite",
    "toolboxes": "toolboxes",
    "mor": "mor",
}
DEFAULT_FULL_PRESET = "default"


def resolve_cmake_preset(target: ImageTarget, *, component: str) -> str:
    normalized = str(component).strip().lower()
    if normalized == "full":
        return target.metadata.get("cmake_full_preset", DEFAULT_FULL_PRESET)
    try:
        default = DEFAULT_COMPONENT_PRESETS[normalized]
    except KeyError as exc:
        known = ", ".join(sorted([*DEFAULT_COMPONENT_PRESETS, "full"]))
        raise ValueError(f"Unsupported CMake preset component '{component}'. Known values: {known}") from exc
    return target.metadata.get(f"cmake_{normalized}_preset", default)
