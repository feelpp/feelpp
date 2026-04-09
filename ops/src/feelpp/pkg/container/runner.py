from __future__ import annotations

from pathlib import Path
import os
import shlex

from ..config import PackagingContext
from ..shell import run
from .args import _containerized_argv, _pass_through_env, publish_uses_aptly, publish_uses_signing
from .constants import (
    CONTAINER_APTLY_CONFIG,
    CONTAINER_APTLY_ROOT,
    CONTAINER_GNUPG_HOME,
    CONTAINER_JOB_ROOT,
    CONTAINER_REPO_ROOT,
)
from .runtime import _bootstrap_apt_packages, default_image_for_context, default_state_root
from .staging import (
    _gpg_bootstrap_commands,
    _host_aptly_config_path,
    _host_gnupg_home,
    _stage_aptly_config,
    _stage_host_feelpp_apt_key,
    _stage_host_gnupg_home,
)


def run_in_docker(
    context: PackagingContext,
    *,
    argv: list[str],
    image: str | None = None,
    state_root: str | Path | None = None,
    dry_run: bool = False,
) -> None:
    image_name = image or default_image_for_context(context)
    host_state_root = Path(state_root).expanduser().resolve() if state_root else default_state_root()
    host_job_root = context.job_root
    host_chroots = host_state_root / "chroots"
    host_aptcache = host_state_root / "aptcache"
    host_aptly = host_state_root / "aptly"
    host_aptly_config = _host_aptly_config_path() if publish_uses_aptly(argv) else None
    host_gnupg_home = _host_gnupg_home() if publish_uses_signing(argv) else None

    for path in (host_state_root, host_chroots, host_aptcache, host_aptly, host_job_root):
        path.mkdir(parents=True, exist_ok=True)
    _stage_host_feelpp_apt_key(job_root=host_job_root)
    if host_aptly_config is not None:
        _stage_aptly_config(host_aptly_config, job_root=host_job_root)
    if host_gnupg_home is not None:
        _stage_host_gnupg_home(host_gnupg_home, job_root=host_job_root)

    command = [
        "docker",
        "run",
        "--rm",
        "--privileged",
        "-v",
        f"{context.repo_root}:{CONTAINER_REPO_ROOT}",
        "-v",
        f"{host_job_root}:{CONTAINER_JOB_ROOT}",
        "-v",
        f"{host_chroots}:/root/pbuilder/chroots",
        "-v",
        f"{host_aptcache}:/var/cache/apt/archives",
        "-v",
        f"{host_aptly}:{CONTAINER_APTLY_ROOT}",
        "-w",
        str(CONTAINER_REPO_ROOT),
    ]

    pass_env = _pass_through_env()
    for key in sorted(pass_env):
        command.extend(["-e", key])
    if host_aptly_config is not None:
        command.extend(["-e", f"FEELPP_APTLY_CONFIG={CONTAINER_APTLY_CONFIG}"])
    command.extend(["-e", "FEELPP_PKG_IN_CONTAINER=1"])

    inner_argv, extra_mounts = _containerized_argv(context, argv)
    for host_path, container_path in extra_mounts:
        command.extend(["-v", f"{host_path}:{container_path}"])
    inner_command = [
        "python3",
        "-m",
        "feelpp.pkg",
        *inner_argv,
    ]
    inner_command_str = " ".join(shlex.quote(token) for token in inner_command)
    bootstrap_lines = [
        "set -e",
        "git config --global --add safe.directory /work || true",
        "git config --global --add safe.directory /work/.git || true",
        "rm -f /etc/apt/apt.conf.d/docker-clean",
        """cat >/etc/apt/apt.conf.d/99feelpp-keep-cache <<'EOF'
APT::Keep-Downloaded-Packages "true";
Binary::apt::APT::Keep-Downloaded-Packages "true";
EOF""",
        f"required_packages={shlex.quote(' '.join(_bootstrap_apt_packages(context)))}",
        'missing_packages=""',
        'for pkg in ${required_packages}; do',
        '  if ! dpkg-query -W -f=\'${Status}\' "${pkg}" 2>/dev/null | grep -q "install ok installed"; then',
        '    missing_packages="${missing_packages} ${pkg}"',
        "  fi",
        "done",
        'if [ -n "${missing_packages}" ]; then',
        "  apt-get update",
        '  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ${missing_packages}',
        "fi",
    ]
    bootstrap_lines.extend(_gpg_bootstrap_commands(host_gnupg_mounted=host_gnupg_home is not None))
    bootstrap_lines.extend(
        [
            f"export PYTHONPATH={shlex.quote(str(CONTAINER_REPO_ROOT / 'ops' / 'src'))}",
            inner_command_str,
        ]
    )
    bootstrap = "\n".join(bootstrap_lines)
    command.extend(
        [
            image_name,
            "bash",
            "-lc",
            bootstrap,
        ]
    )

    run(command, cwd=context.repo_root, env=os.environ.copy(), dry_run=dry_run)
