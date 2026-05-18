from __future__ import annotations

import shlex

from ..config import PackagingContext
from ..shell import run_bash
from .mirrors import pbuilder_mirrorsite, pbuilder_othermirrors


def _run_prepare_base_shell(context: PackagingContext, *, dry_run: bool = False) -> None:
    _run_prepare_base_shell_with_mirror(
        context,
        mirrorsite=pbuilder_mirrorsite(context),
        othermirrors=pbuilder_othermirrors(context),
        dry_run=dry_run,
    )


def _run_prepare_base_shell_with_mirror(
    context: PackagingContext,
    *,
    mirrorsite: str,
    othermirrors: str,
    dry_run: bool = False,
) -> None:
    common_script = context.repo_root / "feelpp" / "tools" / "scripts" / "pkg" / "feelpp_pkg_common.sh"
    command = (
        f"source {shlex.quote(str(common_script))} && "
        f"prepare_feelpp_pbuilder_base {shlex.quote(context.dist)}"
    )
    run_bash(
        command,
        cwd=context.repo_root,
        env=context.shell_env(
            FEELPP_PBUILDER_MIRRORSITE=mirrorsite,
            FEELPP_PBUILDER_OTHERMIRROR=othermirrors,
        ),
        env_overrides={
            "DIST": context.dist,
            "CHANNEL": context.channel,
            "FEELPP_PBUILDER_MIRRORSITE": mirrorsite,
            "FEELPP_PBUILDER_OTHERMIRROR": othermirrors,
        },
        dry_run=dry_run,
    )
