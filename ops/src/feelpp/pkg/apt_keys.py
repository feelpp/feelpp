from __future__ import annotations

from pathlib import Path
import os
import subprocess
from typing import Iterable, Sequence

from .shell import run_checked


CURRENT_APT_SIGNING_KEY = "BD86E2E0A3DA7E56A675D805EF232CA173566681"
LEGACY_APT_SIGNING_KEY = "92C13868485466BB1A6584FA491F361BCEF12211"
DEFAULT_APT_SIGNING_KEY = CURRENT_APT_SIGNING_KEY
DEFAULT_APT_KEYRING_KEYS = (
    LEGACY_APT_SIGNING_KEY,
    CURRENT_APT_SIGNING_KEY,
)


def _dedupe_keys(keys: Iterable[str]) -> list[str]:
    deduped: list[str] = []
    seen: set[str] = set()
    for key in keys:
        normalized = key.strip()
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        deduped.append(normalized)
    return deduped


def _parse_key_list(raw: str) -> list[str]:
    return _dedupe_keys(raw.replace(",", " ").split())


def default_apt_signing_key() -> str:
    return os.getenv("FEELPP_APT_GPG_KEY") or os.getenv("GPG_KEY") or DEFAULT_APT_SIGNING_KEY


def default_apt_keyring_keys(*, primary_key: str | None = None) -> list[str]:
    explicit_keys = os.getenv("FEELPP_APT_GPG_KEYS", "").strip()
    if explicit_keys:
        return _parse_key_list(explicit_keys)

    return _dedupe_keys(
        [
            primary_key or default_apt_signing_key(),
            *DEFAULT_APT_KEYRING_KEYS,
        ]
    )


def export_apt_public_keyring(
    target: Path,
    *,
    key_ids: Sequence[str] | None = None,
) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    resolved_keys = list(key_ids) if key_ids is not None else default_apt_keyring_keys()
    run_checked(
        [
            "gpg",
            "--batch",
            "--yes",
            "--output",
            str(target),
            "--export",
            *resolved_keys,
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
