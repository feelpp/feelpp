from __future__ import annotations

import os
from pathlib import Path
import shlex
import subprocess
from typing import Sequence


def format_command(args: Sequence[str]) -> str:
    return shlex.join([str(arg) for arg in args])


def _print_dry_run_command(
    command: list[str],
    *,
    env_overrides: dict[str, str] | None = None,
) -> None:
    if env_overrides:
        prefix = " ".join(
            f"{key}={shlex.quote(str(value))}" for key, value in env_overrides.items()
        )
        print(f"{prefix} {format_command(command)}")
    else:
        print(format_command(command))


def _effective_env(
    env: dict[str, str] | None,
    env_overrides: dict[str, str] | None,
) -> dict[str, str] | None:
    if not env_overrides:
        return env

    merged = dict(os.environ) if env is None else dict(env)
    merged.update(
        {
            key: str(value)
            for key, value in env_overrides.items()
            if value is not None
        }
    )
    return merged


def run_checked(
    args: Sequence[str],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    env_overrides: dict[str, str] | None = None,
    dry_run: bool = False,
    stdout: int | None = None,
    stderr: int | None = None,
) -> None:
    command = [str(arg) for arg in args]
    if dry_run:
        _print_dry_run_command(command, env_overrides=env_overrides)
        return
    subprocess.run(
        command,
        cwd=cwd,
        env=_effective_env(env, env_overrides),
        check=True,
        stdout=stdout,
        stderr=stderr,
    )


def run_capture(
    args: Sequence[str],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    env_overrides: dict[str, str] | None = None,
    dry_run: bool = False,
    check: bool = True,
    stderr: int | None = None,
) -> str:
    command = [str(arg) for arg in args]
    if dry_run:
        _print_dry_run_command(command, env_overrides=env_overrides)
        return ""
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=_effective_env(env, env_overrides),
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE if stderr is None else stderr,
        text=True,
    )
    if check and completed.returncode != 0:
        raise subprocess.CalledProcessError(
            completed.returncode,
            command,
            output=completed.stdout,
            stderr=completed.stderr if isinstance(completed.stderr, str) else None,
        )
    return completed.stdout


def run_probe(
    args: Sequence[str],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    env_overrides: dict[str, str] | None = None,
    dry_run: bool = False,
    stdout: int | None = None,
    stderr: int | None = None,
) -> bool:
    command = [str(arg) for arg in args]
    if dry_run:
        _print_dry_run_command(command, env_overrides=env_overrides)
        return False
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=_effective_env(env, env_overrides),
        check=False,
        stdout=subprocess.DEVNULL if stdout is None else stdout,
        stderr=subprocess.DEVNULL if stderr is None else stderr,
    )
    return completed.returncode == 0


def run(
    args: Sequence[str],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    env_overrides: dict[str, str] | None = None,
    dry_run: bool = False,
) -> None:
    run_checked(
        args,
        cwd=cwd,
        env=env,
        env_overrides=env_overrides,
        dry_run=dry_run,
    )


def run_bash(
    script: str,
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    env_overrides: dict[str, str] | None = None,
    dry_run: bool = False,
) -> None:
    run_checked(
        ["bash", "-lc", script],
        cwd=cwd,
        env=env,
        env_overrides=env_overrides,
        dry_run=dry_run,
    )
