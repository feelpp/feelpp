from __future__ import annotations

import os

from ..config import PackagingContext


TRUTHY = {"1", "true", "yes", "on"}


def _mirror_site(context: PackagingContext) -> str:
    if context.flavor == "ubuntu":
        return os.getenv(
            "FEELPP_PBUILDER_MIRRORSITE_UBUNTU",
            "http://miroir.univ-lorraine.fr/ubuntu/",
        )
    if context.flavor == "debian":
        return os.getenv(
            "FEELPP_PBUILDER_MIRRORSITE_DEBIAN",
            "http://miroir.univ-lorraine.fr/debian/",
        )
    raise ValueError(f"Unsupported pbuilder mirror flavor: {context.flavor}")


def pbuilder_mirrorsite(context: PackagingContext) -> str:
    return _mirror_site(context)


def _mirror_site_candidates(context: PackagingContext) -> list[str]:
    env_name = f"FEELPP_PBUILDER_MIRRORSITE_{context.flavor.upper()}_FALLBACKS"
    raw = os.getenv(env_name, os.getenv("FEELPP_PBUILDER_MIRRORSITE_FALLBACKS", ""))
    candidates = [_mirror_site(context)]
    for item in raw.replace(",", "|").split("|"):
        value = item.strip()
        if value:
            candidates.append(value)

    deduped: list[str] = []
    seen: set[str] = set()
    for candidate in candidates:
        normalized = candidate.rstrip("/") + "/"
        if normalized in seen:
            continue
        seen.add(normalized)
        deduped.append(normalized)
    return deduped


def local_repo_line(context: PackagingContext) -> str | None:
    packages = context.local_repo_dir / "Packages"
    packages_gz = context.local_repo_dir / "Packages.gz"
    if packages.is_file() or packages_gz.is_file():
        return f"deb [trusted=yes] file://{context.local_repo_dir} ./"
    return None


def pbuilder_othermirrors(
    context: PackagingContext,
    *,
    mirrorsite: str | None = None,
    allow_public_fallback: bool | None = None,
) -> str:
    mirror_root = (mirrorsite or pbuilder_mirrorsite(context)).rstrip("/")
    mirrors: list[str] = []
    local_repo = local_repo_line(context)
    if allow_public_fallback is None:
        allow_public_fallback = (
            os.getenv("FEELPP_PBUILDER_ALLOW_PUBLIC_FEELPP_FALLBACK", "").strip().lower() in TRUTHY
        )
    if local_repo:
        mirrors.append(local_repo)
    elif allow_public_fallback:
        mirrors.append(f"deb http://apt.feelpp.org/{context.flavor}/{context.dist} {context.dist} {context.channel}")

    if context.dist == "focal":
        mirrors.append(f"deb https://apt.kitware.com/{context.flavor}/ {context.dist} main")
    elif context.dist == "jammy":
        mirrors.append(f"deb {mirror_root} {context.dist}-backports main restricted universe multiverse")
    elif context.dist == "bullseye":
        mirrors.append(f"deb {mirror_root} {context.dist}-backports main")

    return "|".join(mirrors)
