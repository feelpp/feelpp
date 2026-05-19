from __future__ import annotations

from collections.abc import Sequence


def normalize_platform_overrides(values: Sequence[str] | None) -> list[str] | None:
    if not values:
        return None

    normalized: list[str] = []
    seen: set[str] = set()
    for raw_value in values:
        for item in str(raw_value).replace("\n", ",").split(","):
            platform = item.strip()
            if not platform or platform in seen:
                continue
            seen.add(platform)
            normalized.append(platform)
    return normalized or None


def resolve_target_platforms(
    default_platforms: Sequence[str],
    overrides: Sequence[str] | None = None,
) -> list[str]:
    normalized = normalize_platform_overrides(overrides)
    if normalized is not None:
        return normalized
    return list(default_platforms)
