from __future__ import annotations

from pathlib import Path


CONTAINER_REPO_ROOT = Path("/work")
CONTAINER_JOB_ROOT = Path("/tmp/feelpp-pkg-job")
CONTAINER_EXTRA_ROOT = Path("/tmp/feelpp-pkg-extra")
CONTAINER_PBUILDER_ROOT = Path("/root/pbuilder")
CONTAINER_APTLY_ROOT = Path("/srv/aptly")
CONTAINER_APTLY_CONFIG = CONTAINER_JOB_ROOT / "aptly.conf"
CONTAINER_GNUPG_HOME = Path("/root/.gnupg")
CONTAINER_STAGED_GNUPG_HOME = CONTAINER_JOB_ROOT / "gnupg-home"

DEFAULT_DOCKER_IMAGE = "ghcr.io/feelpp/pkg-env:ubuntu-24.04"
DEFAULT_DOCKER_IMAGES_BY_CONTEXT: dict[tuple[str, str], str] = {
    ("ubuntu", "noble"): "ghcr.io/feelpp/pkg-env:ubuntu-24.04",
    ("ubuntu", "resolute"): "ghcr.io/feelpp/pkg-env:ubuntu-26.04",
    ("debian", "trixie"): "ghcr.io/feelpp/pkg-env:debian-trixie",
}

DEFAULT_APT_SIGNING_KEY = "BD86E2E0A3DA7E56A675D805EF232CA173566681"
DEFAULT_STAGED_APT_KEYRING = "feelpp-archive-keyring.gpg"
