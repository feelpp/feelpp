from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import os


UBUNTU_DISTS = {"focal", "jammy", "lunar", "mantic", "noble", "resolute"}
DEBIAN_DISTS = {"bullseye", "bookworm", "trixie", "testing", "sid"}
FEDORA_DISTS = {"fedora-42"}


def discover_repo_root() -> Path:
    env_root = os.getenv("FEELPP_REPO_ROOT")
    if env_root:
        return Path(env_root).expanduser().resolve()

    current = Path.cwd().resolve()
    for candidate in (current, *current.parents):
        if (candidate / "packaging" / "debian").is_dir() and (
            candidate / "feelpp" / "tools" / "scripts" / "pkg" / "feelpp_pkg.sh"
        ).is_file():
            return candidate

    return Path(__file__).resolve().parents[4]


def detect_branch(explicit_branch: str | None = None) -> str:
    return (
        explicit_branch
        or os.getenv("BRANCH")
        or os.getenv("BUILDKITE_BRANCH")
        or os.getenv("GITHUB_REF_NAME")
        or "develop"
    )


def detect_channel(branch: str, explicit_channel: str | None = None) -> str:
    if explicit_channel:
        return explicit_channel
    if branch in {"main", "master"}:
        return "stable"
    return "latest"


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


def detect_job_id(explicit_job_id: str | None = None) -> str:
    if explicit_job_id:
        return explicit_job_id
    if os.getenv("FEELPP_PKG_JOB_ID"):
        return os.environ["FEELPP_PKG_JOB_ID"]
    run_id = os.getenv("GITHUB_RUN_ID", "local")
    run_attempt = os.getenv("GITHUB_RUN_ATTEMPT", "0")
    job = os.getenv("GITHUB_JOB", "pkg")
    return f"{run_id}-{run_attempt}-{job}"


@dataclass(frozen=True)
class PackagingContext:
    repo_root: Path
    dist: str
    flavor: str
    branch: str
    channel: str
    job_id: str
    job_root: Path
    local_repo_dir: Path
    artifacts_dir: Path
    results_dir: Path
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
    ) -> "PackagingContext":
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
            dist=resolved_dist,
            flavor=resolved_flavor,
            branch=resolved_branch,
            channel=resolved_channel,
            job_id=resolved_job_id,
            job_root=resolved_job_root,
            local_repo_dir=local_repo_dir,
            artifacts_dir=artifacts_dir,
            results_dir=results_dir,
            pbuilder_root=pbuilder_root,
            pbuilder_config=pbuilder_config,
            pbuilder_source_hookdir=pbuilder_source_hookdir,
            pbuilder_runtime_hookdir=pbuilder_runtime_hookdir,
            pbuilder_keyrings_dir=pbuilder_keyrings_dir,
        )

    @property
    def manifest_path(self) -> Path:
        return self.repo_root / "packaging" / "manifest" / "components.toml"

    @property
    def job_manifest_path(self) -> Path:
        return self.job_root / "feelpp-pkg-job.json"

    def component_results_dir(self, component: str) -> Path:
        return self.results_dir / component

    def shell_env(self, **extra: str) -> dict[str, str]:
        env = dict(os.environ)
        env.update(
            {
                "FEELPP_REPO_ROOT": str(self.repo_root),
                "DIST": self.dist,
                "FLAVOR": self.flavor,
                "BRANCH": self.branch,
                "CHANNEL": self.channel,
                "FEELPP_PKG_JOB_ID": self.job_id,
                "FEELPP_PKG_JOB_ROOT": str(self.job_root),
                "FEELPP_PKG_LOCAL_REPO_DIR": str(self.local_repo_dir),
                "FEELPP_PKG_ARTIFACTS_DIR": str(self.artifacts_dir),
                "FEELPP_PKG_RESULTS_DIR": str(self.results_dir),
                "FEELPP_PBUILDER_ROOT": str(self.pbuilder_root),
                "FEELPP_PBUILDER_CONFIG": str(self.pbuilder_config),
                "FEELPP_PBUILDER_SOURCE_HOOKDIR": str(self.pbuilder_source_hookdir),
                "FEELPP_PBUILDER_HOOKDIR": str(self.pbuilder_runtime_hookdir),
                "FEELPP_PBUILDER_KEYRINGS_DIR": str(self.pbuilder_keyrings_dir),
            }
        )
        env.update({key: value for key, value in extra.items() if value is not None})
        return env

    def as_dict(self) -> dict[str, str]:
        data = asdict(self)
        return {key: str(value) for key, value in data.items()}
