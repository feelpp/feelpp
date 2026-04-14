from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import os
import tempfile

from ..config import PackagingContext
from ..shell import run_capture, run_checked, run_probe
from ..workspace import read_job_manifest


@dataclass(frozen=True)
class RepoPackageRef:
    raw: str
    name: str
    version: str
    architecture: str


def assert_publish_ready(context: PackagingContext) -> None:
    job_manifest = read_job_manifest(context)
    if not job_manifest:
        return

    plan = job_manifest.get("plan")
    if not isinstance(plan, dict):
        return

    components = plan.get("components")
    if not isinstance(components, list):
        return

    expected = [
        component.get("name")
        for component in components
        if isinstance(component, dict) and isinstance(component.get("name"), str)
    ]
    if not expected:
        return

    built_components = {
        component
        for component in job_manifest.get("built_components", [])
        if isinstance(component, str)
    }
    missing = [component for component in expected if component not in built_components]
    if not missing:
        return

    raise RuntimeError(
        "Refusing publish because the recorded build chain is incomplete. "
        f"Missing built components: {', '.join(missing)}."
    )


def snapshot_id() -> str:
    explicit = os.getenv("FEELPP_PKG_SNAPSHOT_ID")
    if explicit:
        return explicit
    run_id = os.getenv("GITHUB_RUN_ID", "local")
    run_attempt = os.getenv("GITHUB_RUN_ATTEMPT", "0")
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")
    return f"{run_id}-{run_attempt}-{timestamp}"


def aptly_base_command() -> list[str]:
    command = ["aptly"]
    aptly_config = os.getenv("FEELPP_APTLY_CONFIG")
    if aptly_config:
        command.append(f"-config={aptly_config}")
    return command


def has_binary_packages(input_dir: Path) -> bool:
    return bool(binary_package_paths(input_dir))


def binary_package_paths(input_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in input_dir.rglob("*")
        if path.is_file() and path.suffix in {".deb", ".udeb"}
    )


def build_publish_args(*, passphrase_file: Path | None) -> list[str]:
    publish_args = ["-force-overwrite"]
    if os.getenv("FEELPP_APTLY_SKIP_SIGNING", "").lower() == "true":
        publish_args.append("-skip-signing")
        return publish_args

    if passphrase_file is not None:
        publish_args.extend(["-batch", f"-passphrase-file={passphrase_file}"])
    gpg_key = os.getenv("GPG_KEY")
    if gpg_key:
        publish_args.append(f"-gpg-key={gpg_key}")
    return publish_args


def create_passphrase_file(context: PackagingContext) -> Path | None:
    passphrase = os.getenv("GPG_PASSPHRASE")
    if not passphrase or os.getenv("FEELPP_APTLY_SKIP_SIGNING", "").lower() == "true":
        return None
    context.job_root.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=context.job_root,
        prefix="aptly-passphrase.",
        delete=False,
    ) as handle:
        handle.write(passphrase)
        path = Path(handle.name)
    path.chmod(0o600)
    return path


def publish_identity(context: PackagingContext) -> tuple[str, str, str]:
    publish_distribution = os.getenv("FEELPP_APTLY_PUBLISH_DISTRIBUTION", context.dist)
    publish_component = os.getenv("FEELPP_APTLY_PUBLISH_COMPONENT", context.channel)
    publish_target = os.getenv(
        "FEELPP_APTLY_PUBLISH_TARGET",
        f"s3:apt.feelpp.org:{context.flavor}/{context.dist}",
    )
    return publish_distribution, publish_component, publish_target


def repo_name(context: PackagingContext) -> str:
    return f"feelpp-{context.dist}-{context.channel}"


def publish_snapshot_from_repo(
    context: PackagingContext,
    *,
    env: dict[str, str],
    snapshot_name: str,
    publish_distribution: str,
    publish_component: str,
    publish_target: str,
    publish_args: list[str],
    dry_run: bool,
) -> None:
    aptly = aptly_base_command()
    run_checked(
        aptly + ["snapshot", "create", snapshot_name, "from", "repo", repo_name(context)],
        cwd=context.repo_root,
        env=env,
        dry_run=dry_run,
    )

    publish_exists = False if dry_run else run_probe(
        aptly + ["publish", "show", publish_distribution, publish_target],
        cwd=context.repo_root,
        env=env,
    )

    if publish_exists:
        run_checked(
            aptly
            + [
                "publish",
                "switch",
                *publish_args,
                f"-component={publish_component}",
                publish_distribution,
                publish_target,
                snapshot_name,
            ],
            cwd=context.repo_root,
            env=env,
            dry_run=dry_run,
        )
    else:
        run_checked(
            aptly
            + [
                "publish",
                "snapshot",
                *publish_args,
                f"-distribution={publish_distribution}",
                f"-component={publish_component}",
                snapshot_name,
                publish_target,
            ],
            cwd=context.repo_root,
            env=env,
            dry_run=dry_run,
        )
