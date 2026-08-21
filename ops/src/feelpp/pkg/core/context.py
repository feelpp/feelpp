from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import os


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

    return Path(__file__).resolve().parents[5]


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
class WorkspaceContext:
    repo_root: Path
    branch: str
    channel: str
    job_id: str
    job_root: Path
    artifacts_dir: Path
    results_dir: Path

    @classmethod
    def create(
        cls,
        *,
        repo_root: str | Path | None = None,
        branch: str | None = None,
        channel: str | None = None,
        job_id: str | None = None,
        job_root: str | Path | None = None,
    ) -> "WorkspaceContext":
        resolved_repo_root = Path(repo_root).expanduser().resolve() if repo_root else discover_repo_root()
        resolved_branch = detect_branch(branch)
        resolved_channel = detect_channel(resolved_branch, channel or os.getenv("CHANNEL"))
        resolved_job_id = detect_job_id(job_id)

        default_job_root = resolved_repo_root / ".cache" / "feelpp-pkg" / "jobs" / resolved_job_id
        resolved_job_root = Path(
            os.getenv("FEELPP_PKG_JOB_ROOT", str(job_root or default_job_root))
        ).expanduser().resolve()

        artifacts_dir = Path(
            os.getenv("FEELPP_PKG_ARTIFACTS_DIR", str(resolved_job_root / "artifacts"))
        ).expanduser().resolve()
        results_dir = Path(
            os.getenv("FEELPP_PKG_RESULTS_DIR", str(resolved_job_root / "results"))
        ).expanduser().resolve()

        return cls(
            repo_root=resolved_repo_root,
            branch=resolved_branch,
            channel=resolved_channel,
            job_id=resolved_job_id,
            job_root=resolved_job_root,
            artifacts_dir=artifacts_dir,
            results_dir=results_dir,
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
                "BRANCH": self.branch,
                "CHANNEL": self.channel,
                "FEELPP_PKG_JOB_ID": self.job_id,
                "FEELPP_PKG_JOB_ROOT": str(self.job_root),
                "FEELPP_PKG_ARTIFACTS_DIR": str(self.artifacts_dir),
                "FEELPP_PKG_RESULTS_DIR": str(self.results_dir),
            }
        )
        env.update({key: value for key, value in extra.items() if value is not None})
        return env

    def as_dict(self) -> dict[str, str]:
        data = asdict(self)
        return {key: str(value) for key, value in data.items()}
