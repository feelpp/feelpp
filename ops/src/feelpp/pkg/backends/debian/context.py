from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os

from ...core.context import WorkspaceContext, detect_branch, detect_channel, detect_job_id, discover_repo_root


UBUNTU_DISTS = {"focal", "jammy", "lunar", "mantic", "noble", "resolute"}
DEBIAN_DISTS = {"bullseye", "bookworm", "trixie", "testing", "sid"}
FEDORA_DISTS = {"fedora-42"}


def detect_flavor(dist: str, explicit_flavor: str | None = None) -> str:
    if explicit_flavor:
        return explicit_flavor
    if dist in UBUNTU_DISTS:
        return "ubuntu"
    if dist in DEBIAN_DISTS:
        return "debian"
    if dist in FEDORA_DISTS:
        return "fedora"
    raise ValueError(f"Unsupported DIST: {dist}")


@dataclass(frozen=True)
class DebianPackagingContext(WorkspaceContext):
    dist: str
    flavor: str
    local_repo_dir: Path
    pbuilder_root: Path
    pbuilder_config: Path
    pbuilder_source_hookdir: Path
    pbuilder_runtime_hookdir: Path
    pbuilder_keyrings_dir: Path

    @classmethod
    def create(
        cls,
        *,
        repo_root: str | Path | None = None,
        dist: str | None = None,
        flavor: str | None = None,
        branch: str | None = None,
        channel: str | None = None,
        job_id: str | None = None,
        job_root: str | Path | None = None,
    ) -> "DebianPackagingContext":
        resolved_repo_root = Path(repo_root).expanduser().resolve() if repo_root else discover_repo_root()
        resolved_dist = dist or os.getenv("DIST") or "noble"
        resolved_branch = detect_branch(branch)
        resolved_channel = detect_channel(resolved_branch, channel or os.getenv("CHANNEL"))
        resolved_flavor = detect_flavor(resolved_dist, flavor or os.getenv("FLAVOR"))
        resolved_job_id = detect_job_id(job_id)

        default_job_root = (
            Path.home() / "pbuilder" / "jobs" / resolved_job_id / resolved_flavor / resolved_dist
        )
        resolved_job_root = Path(
            os.getenv("FEELPP_PKG_JOB_ROOT", str(job_root or default_job_root))
        ).expanduser().resolve()

        local_repo_dir = Path(
            os.getenv("FEELPP_PKG_LOCAL_REPO_DIR", str(resolved_job_root / "local-repo"))
        ).expanduser().resolve()
        artifacts_dir = Path(
            os.getenv("FEELPP_PKG_ARTIFACTS_DIR", str(resolved_job_root / "artifacts"))
        ).expanduser().resolve()
        results_dir = Path(
            os.getenv("FEELPP_PKG_RESULTS_DIR", str(resolved_job_root / "results"))
        ).expanduser().resolve()

        pbuilder_root = Path(
            os.getenv(
                "FEELPP_PBUILDER_ROOT",
                str(
                    Path.home()
                    / "pbuilder"
                    / "chroots"
                    / resolved_flavor
                    / resolved_dist
                    / resolved_channel
                ),
            )
        ).expanduser().resolve()
        pbuilder_dir = Path(
            os.getenv(
                "FEELPP_PBUILDER_DIR",
                str(resolved_repo_root / "packaging" / "pbuilder"),
            )
        ).expanduser().resolve()
        pbuilder_config = Path(
            os.getenv("FEELPP_PBUILDER_CONFIG", str(pbuilder_dir / "pbuilderrc"))
        ).expanduser().resolve()
        pbuilder_source_hookdir = Path(
            os.getenv("FEELPP_PBUILDER_SOURCE_HOOKDIR", str(pbuilder_dir / "hooks"))
        ).expanduser().resolve()
        pbuilder_runtime_hookdir = Path(
            os.getenv("FEELPP_PBUILDER_HOOKDIR", str(resolved_job_root / "pbuilder" / "hooks"))
        ).expanduser().resolve()
        pbuilder_keyrings_dir = Path(
            os.getenv("FEELPP_PBUILDER_KEYRINGS_DIR", str(resolved_job_root / "pbuilder" / "keyrings"))
        ).expanduser().resolve()

        return cls(
            repo_root=resolved_repo_root,
            branch=resolved_branch,
            channel=resolved_channel,
            job_id=resolved_job_id,
            job_root=resolved_job_root,
            artifacts_dir=artifacts_dir,
            results_dir=results_dir,
            dist=resolved_dist,
            flavor=resolved_flavor,
            local_repo_dir=local_repo_dir,
            pbuilder_root=pbuilder_root,
            pbuilder_config=pbuilder_config,
            pbuilder_source_hookdir=pbuilder_source_hookdir,
            pbuilder_runtime_hookdir=pbuilder_runtime_hookdir,
            pbuilder_keyrings_dir=pbuilder_keyrings_dir,
        )

    def shell_env(self, **extra: str) -> dict[str, str]:
        env = super().shell_env(
            DIST=self.dist,
            FLAVOR=self.flavor,
            FEELPP_PKG_LOCAL_REPO_DIR=str(self.local_repo_dir),
            FEELPP_PBUILDER_ROOT=str(self.pbuilder_root),
            FEELPP_PBUILDER_CONFIG=str(self.pbuilder_config),
            FEELPP_PBUILDER_SOURCE_HOOKDIR=str(self.pbuilder_source_hookdir),
            FEELPP_PBUILDER_HOOKDIR=str(self.pbuilder_runtime_hookdir),
            FEELPP_PBUILDER_KEYRINGS_DIR=str(self.pbuilder_keyrings_dir),
        )
        env.update({key: value for key, value in extra.items() if value is not None})
        return env
