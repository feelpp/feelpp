from __future__ import annotations

import os
from pathlib import Path
import subprocess

import pytest
import yaml


WORKFLOW = Path(__file__).resolve().parents[2] / ".github/workflows/ci-ops.yml"


def _git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def _commit(repo: Path) -> str:
    _git(repo, "add", ".")
    _git(
        repo, "-c", "user.name=CI Test", "-c", "user.email=ci@example.invalid",
        "-c", "commit.gpgsign=false", "commit", "-qm", "Test change",
    )
    return _git(repo, "rev-parse", "HEAD")


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    _git(tmp_path, "init", "-q")
    for name in (
        "ops/README.md",
        "ops/src/tool.py",
        "packaging/debian/feelpp/noble/debian/changelog",
    ):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("Initial content\n", encoding="utf-8")
    _commit(tmp_path)
    return tmp_path


def _classify(repo: Path, base: str) -> str:
    workflow = yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))
    step = next(
        step for step in workflow["jobs"]["classify"]["steps"]
        if step.get("id") == "classify"
    )
    output = repo / "github-output"
    env = {
        **os.environ,
        "BASE_SHA": base,
        "HEAD_SHA": _git(repo, "rev-parse", "HEAD"),
        "GITHUB_OUTPUT": str(output),
    }
    subprocess.run(["bash", "-e", "-o", "pipefail", "-c", step["run"]],
                   cwd=repo, env=env, check=True, capture_output=True, text=True)
    return output.read_text(encoding="utf-8")


def test_metadata_diff_does_not_run_structural_checks(repo: Path) -> None:
    base = _git(repo, "rev-parse", "HEAD")
    changelog = "packaging/debian/feelpp/noble/debian/changelog"
    (repo / changelog).write_text("Updated metadata\n", encoding="utf-8")
    _commit(repo)

    output = _classify(repo, base)

    assert "structural=false\n" in output
    assert f"changed_changelogs<<EOF\n{changelog}\nEOF\n" in output


def test_structural_diff_runs_ops_checks(repo: Path) -> None:
    base = _git(repo, "rev-parse", "HEAD")
    (repo / "ops/src/tool.py").write_text("Updated tool\n", encoding="utf-8")
    _commit(repo)

    output = _classify(repo, base)

    assert "structural=true\n" in output
    assert "changed_changelogs<<EOF\nEOF\n" in output


@pytest.mark.parametrize("base", ["1" * 40, "0" * 40, ""])
def test_unavailable_base_validates_current_tree(repo: Path, base: str) -> None:
    output = _classify(repo, base)

    assert "structural=true\n" in output
    assert "packaging/debian/feelpp/noble/debian/changelog\n" in output
