from __future__ import annotations

from pathlib import Path

from ..config import PackagingContext
from ..pbuilder import pbuilder_mirrorsite, pbuilder_othermirrors
from ..shell import run


def run_pbuilder_build(
    context: PackagingContext,
    dsc_path: Path,
    result_dir: Path,
    *,
    allow_public_fallback: bool = False,
    dry_run: bool = False,
) -> None:
    result_dir.mkdir(parents=True, exist_ok=True)
    buildplace = context.pbuilder_root / "build"
    buildplace.mkdir(parents=True, exist_ok=True)
    auth_script = context.repo_root / "feelpp" / "tools" / "scripts" / "pkg" / "feelpp_pkg_sudo_auth.sh"

    mirrorsite = pbuilder_mirrorsite(context)
    othermirrors = pbuilder_othermirrors(
        context,
        mirrorsite=mirrorsite,
        allow_public_fallback=allow_public_fallback,
    )
    env_overrides = {
        "CHANNEL": context.channel,
        "DIST": context.dist,
        "FEELPP_PBUILDER_ALLOW_PUBLIC_FEELPP_FALLBACK": "true" if allow_public_fallback else "false",
        "FEELPP_PBUILDER_BINDMOUNTS": str(context.local_repo_dir),
        "FLAVOR": context.flavor,
        "MIRRORSITE": mirrorsite,
        "OTHERMIRROR": othermirrors,
        "PBUILDAUTH": str(auth_script),
        "PBUILDFOLDER": str(context.pbuilder_root),
    }
    command = [
        "pbuilder-dist",
        context.dist,
        "build",
        "--configfile",
        str(context.pbuilder_config),
        "--buildresult",
        str(result_dir),
        "--buildplace",
        str(buildplace),
        "--hookdir",
        str(context.pbuilder_runtime_hookdir),
        "--mirror",
        mirrorsite,
        "--othermirror",
        othermirrors,
        str(dsc_path),
    ]
    run(
        command,
        cwd=context.repo_root,
        env=context.shell_env(),
        env_overrides=env_overrides,
        dry_run=dry_run,
    )
