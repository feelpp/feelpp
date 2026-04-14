from __future__ import annotations

from pathlib import Path
import os

from ..config import PackagingContext
from .common import (
    aptly_base_command,
    assert_publish_ready,
    binary_package_paths,
    build_publish_args,
    create_passphrase_file,
    publish_identity,
    publish_snapshot_from_repo,
    repo_name,
    run_checked,
    run_probe,
    snapshot_id,
)


def _dry_run_publish(
    *,
    cwd: Path,
    env: dict[str, str],
    repo_name_value: str,
    snapshot_name: str,
    publish_distribution: str,
    publish_component: str,
    publish_target: str,
    input_dir: Path,
    binary_packages: list[Path],
    has_binary_packages_value: bool,
    publish_args: list[str],
) -> None:
    base = aptly_base_command()
    run_checked(
        base
        + [
            "repo",
            "create",
            f"-distribution={publish_distribution}",
            f"-component={publish_component}",
            repo_name_value,
        ],
        cwd=cwd,
        env=env,
        dry_run=True,
    )
    if has_binary_packages_value:
        run_checked(
            base
            + ["repo", "add", "-force-replace", repo_name_value]
            + [str(path) for path in binary_packages],
            cwd=cwd,
            env=env,
            dry_run=True,
        )
    else:
        print(f"No binary packages found in {input_dir}, skipping aptly repo add")
    run_checked(
        base + ["snapshot", "create", snapshot_name, "from", "repo", repo_name_value],
        cwd=cwd,
        env=env,
        dry_run=True,
    )
    run_checked(
        base
        + [
            "publish",
            "snapshot",
            *publish_args,
            f"-distribution={publish_distribution}",
            f"-component={publish_component}",
            snapshot_name,
            publish_target,
        ],
        cwd=cwd,
        env=env,
        dry_run=True,
    )
    print(f"Published snapshot: {snapshot_name}")


def publish_snapshot(
    context: PackagingContext,
    *,
    input_dir: Path | None = None,
    dry_run: bool = False,
) -> None:
    assert_publish_ready(context)

    effective_input_dir = (input_dir or context.artifacts_dir).expanduser().resolve()
    env = context.shell_env(FEELPP_PKG_PUBLISH_INPUT_DIR=str(effective_input_dir))

    repo_name_value = repo_name(context)
    snapshot_name = f"{repo_name_value}-snapshot-{snapshot_id()}"
    publish_distribution, publish_component, publish_target = publish_identity(context)

    if dry_run:
        binary_packages = binary_package_paths(effective_input_dir) if effective_input_dir.is_dir() else []
        publish_args = build_publish_args(
            passphrase_file=Path("/tmp/aptly-passphrase") if os.getenv("GPG_PASSPHRASE") else None
        )
        _dry_run_publish(
            cwd=context.repo_root,
            env=env,
            repo_name_value=repo_name_value,
            snapshot_name=snapshot_name,
            publish_distribution=publish_distribution,
            publish_component=publish_component,
            publish_target=publish_target,
            input_dir=effective_input_dir,
            binary_packages=binary_packages,
            has_binary_packages_value=bool(binary_packages),
            publish_args=publish_args,
        )
        return

    if not effective_input_dir.is_dir():
        raise FileNotFoundError(f"Publish input directory not found: {effective_input_dir}")

    passphrase_file = create_passphrase_file(context)
    try:
        publish_args = build_publish_args(passphrase_file=passphrase_file)
        aptly = aptly_base_command()
        binary_packages = binary_package_paths(effective_input_dir)

        repo_exists = run_probe(
            aptly + ["repo", "show", repo_name_value],
            cwd=context.repo_root,
            env=env,
        )
        if not repo_exists:
            run_checked(
                aptly
                + [
                    "repo",
                    "create",
                    f"-distribution={publish_distribution}",
                    f"-component={publish_component}",
                    repo_name_value,
                ],
                cwd=context.repo_root,
                env=env,
                dry_run=False,
            )

        if binary_packages:
            run_checked(
                aptly
                + ["repo", "add", "-force-replace", repo_name_value]
                + [str(path) for path in binary_packages],
                cwd=context.repo_root,
                env=env,
                dry_run=False,
            )
        else:
            print(f"No binary packages found in {effective_input_dir}, skipping aptly repo add")
        publish_snapshot_from_repo(
            context,
            env=env,
            snapshot_name=snapshot_name,
            publish_distribution=publish_distribution,
            publish_component=publish_component,
            publish_target=publish_target,
            publish_args=publish_args,
            dry_run=False,
        )

        print(f"Published snapshot: {snapshot_name}")
    finally:
        if passphrase_file is not None and passphrase_file.exists():
            passphrase_file.unlink()
