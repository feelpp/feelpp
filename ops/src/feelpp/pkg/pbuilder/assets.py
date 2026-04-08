from __future__ import annotations

import base64
from pathlib import Path
import os
import shutil
import shlex
import subprocess

from ..config import PackagingContext
from ..shell import run_checked
from ..workspace import ensure_workspace
from .constants import STAGED_APT_KEYRING_NAME
from .seed import _batched, _builddeps_batch_size, _load_builddeps


DEFAULT_APT_ACQUIRE_RETRIES = 8
DEFAULT_APT_HTTP_TIMEOUT = 30
DEFAULT_APT_SIGNING_KEY = "BD86E2E0A3DA7E56A675D805EF232CA173566681"


def resolve_packaging_tree(context: PackagingContext, component: str) -> tuple[Path, Path]:
    component_dir = context.repo_root / "packaging" / "debian" / component
    dist_dir = component_dir / context.dist

    if not component_dir.is_dir():
        raise FileNotFoundError(
            f"Packaging metadata not found for component {component}: {component_dir}"
        )
    if not dist_dir.is_dir():
        raise FileNotFoundError(
            f"Packaging metadata not found for {component}/{context.dist}: {dist_dir}"
        )
    return component_dir, dist_dir


def _generate_builddeps_hook(context: PackagingContext) -> str:
    builddeps = _load_builddeps(context)
    acquire_retries = DEFAULT_APT_ACQUIRE_RETRIES
    http_timeout = DEFAULT_APT_HTTP_TIMEOUT
    batch_commands = "\n".join(
        f"install_batch {' '.join(shlex.quote(entry) for entry in batch)}"
        for batch in _batched(builddeps, _builddeps_batch_size())
    )
    return f"""#!/bin/sh

set -eu

apt_update() {{
    apt-get -o Acquire::Retries={acquire_retries} -o Acquire::http::Timeout={http_timeout} update -yq
}}

install_batch() {{
    attempt=1
    max_attempts=5
    while true; do
        if apt-get -o Acquire::Retries={acquire_retries} -o Acquire::http::Timeout={http_timeout} -o APT::Install-Recommends=false -o APT::Install-Suggests=false install -yq --no-install-recommends --fix-missing "$@"; then
            return 0
        fi
        if [ "${{attempt}}" -ge "${{max_attempts}}" ]; then
            return 1
        fi
        attempt=$((attempt + 1))
        apt_update
    done
}}

apt_update
{batch_commands}

case "${{DISTRIBUTION}}" in
    bullseye)
        attempt=1
        max_attempts=5
        while true; do
            if apt-get -o Acquire::Retries={acquire_retries} -o Acquire::http::Timeout={http_timeout} install -yq --fix-missing -t bullseye-backports cmake; then
                break
            fi
            if [ "${{attempt}}" -ge "${{max_attempts}}" ]; then
                exit 1
            fi
            attempt=$((attempt + 1))
            apt_update
        done
        ;;
esac
"""


def _feelpp_apt_signing_key() -> str:
    return os.getenv("FEELPP_APT_GPG_KEY") or os.getenv("GPG_KEY") or DEFAULT_APT_SIGNING_KEY


def refresh_feelpp_keyring(target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        run_checked(
            [
                "gpg",
                "--batch",
                "--yes",
                "--output",
                str(target),
                "--export",
                _feelpp_apt_signing_key(),
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return
    if target.is_file():
        target.chmod(0o644)


def _generate_keyrings_hook(context: PackagingContext) -> str:
    lines = [
        "#!/bin/sh",
        "",
        "set -eu",
        "",
        "install -d -m 0755 /etc/apt/trusted.gpg.d",
    ]

    for keyring in sorted(context.pbuilder_keyrings_dir.glob("*.gpg")):
        payload = base64.b64encode(keyring.read_bytes()).decode("ascii")
        target_path = f"/etc/apt/trusted.gpg.d/{keyring.name}"
        lines.extend(
            [
                f"base64 -d >{shlex.quote(target_path)} <<'EOF_{keyring.stem.upper()}'",
                payload,
                f"EOF_{keyring.stem.upper()}",
                f"chmod 0644 {shlex.quote(target_path)}",
            ]
        )

    lines.append("")
    return "\n".join(lines)


def prepare_runtime_assets(context: PackagingContext) -> None:
    ensure_workspace(context)
    context.pbuilder_keyrings_dir.mkdir(parents=True, exist_ok=True)
    context.pbuilder_runtime_hookdir.mkdir(parents=True, exist_ok=True)

    if not context.pbuilder_source_hookdir.is_dir():
        raise FileNotFoundError(
            f"pbuilder source hook directory not found: {context.pbuilder_source_hookdir}"
        )

    keyring_source_dir = context.pbuilder_source_hookdir / "keyrings"
    if not keyring_source_dir.is_dir():
        raise FileNotFoundError(f"pbuilder keyring directory not found: {keyring_source_dir}")

    for existing in context.pbuilder_keyrings_dir.glob("*.gpg"):
        existing.unlink()

    for encoded in sorted(keyring_source_dir.glob("*.gpg.b64")):
        keyring_name = encoded.name.removesuffix(".b64")
        decoded = base64.b64decode(encoded.read_text(encoding="utf-8"))
        target = context.pbuilder_keyrings_dir / keyring_name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(decoded)
        target.chmod(0o644)

    runtime_feelpp_keyring = context.pbuilder_keyrings_dir / "feelpp.gpg"
    staged_keyring = context.job_root / STAGED_APT_KEYRING_NAME
    if staged_keyring.is_file():
        shutil.copy2(staged_keyring, runtime_feelpp_keyring)
        runtime_feelpp_keyring.chmod(0o644)
    else:
        refresh_feelpp_keyring(runtime_feelpp_keyring)

    for existing in context.pbuilder_runtime_hookdir.iterdir():
        if existing.is_file() or existing.is_symlink():
            existing.unlink()
        elif existing.is_dir():
            shutil.rmtree(existing)

    for hook in sorted(context.pbuilder_source_hookdir.glob("[A-Z0-9]*-*")):
        if hook.name in {"E10-feelpp-builddeps", "G10-feelpp-keyrings"}:
            continue
        target = context.pbuilder_runtime_hookdir / hook.name
        shutil.copy2(hook, target)
        target.chmod(0o755)

    keyrings_hook = context.pbuilder_runtime_hookdir / "G10-feelpp-keyrings"
    keyrings_hook.write_text(_generate_keyrings_hook(context), encoding="utf-8")
    keyrings_hook.chmod(0o755)

    builddeps_hook = context.pbuilder_runtime_hookdir / "E10-feelpp-builddeps"
    builddeps_hook.write_text(_generate_builddeps_hook(context), encoding="utf-8")
    builddeps_hook.chmod(0o755)
