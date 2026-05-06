from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterable

from feelpp.pkg.core.context import discover_repo_root
from feelpp.pkg.oci.catalog import ImageTarget, get_image_target, load_image_profile_config
from feelpp.pkg.oci.common import oci_image_ref


DEFAULT_IMAGE_PROFILE = "images"
DEFAULT_DEVCONTAINER_GROUP = "env-current"
DEFAULT_DEVCONTAINER_TARGET = "debian:trixie"
DEFAULT_REMOTE_USER = "feelpp"
DEFAULT_BUILD_JOBS = 25
DEVCONTAINER_DIR = ".devcontainer"
ROOT_DEVCONTAINER_PATH = Path(DEVCONTAINER_DIR) / "devcontainer.json"
CMAKE_USER_PRESETS_PATH = Path("CMakeUserPresets.json")

VSCODE_EXTENSIONS = [
    "ms-vscode.cpptools",
    "ms-vscode.cmake-tools",
    "josetr.cmake-language-support-vscode",
    "asciidoctor.asciidoctor-vscode",
    "ms-python.python",
    "ms-toolsai.jupyter",
    "ms-vsliveshare.vsliveshare",
    "shardulm94.trailing-spaces",
    "redhat.vscode-yaml",
]


@dataclass(frozen=True)
class DevcontainerProfile:
    id: str
    target: str
    name: str
    image: str
    cmake_preset: str
    inherits_preset: str
    distribution_env: str
    build_dir: str
    install_dir: str
    remote_user: str = DEFAULT_REMOTE_USER
    build_jobs: int = DEFAULT_BUILD_JOBS

    @property
    def relative_path(self) -> Path:
        return Path(DEVCONTAINER_DIR) / self.id / "devcontainer.json"

    def as_row(self) -> dict[str, str]:
        return {
            "id": self.id,
            "target": self.target,
            "image": self.image,
            "cmake_preset": self.cmake_preset,
            "inherits_preset": self.inherits_preset,
            "path": str(self.relative_path),
        }


def _repo_root(raw: str | Path | None) -> Path:
    return Path(raw).expanduser().resolve() if raw else discover_repo_root()


def _json(payload: object) -> str:
    return json.dumps(payload, indent=3) + "\n"


def _title(value: str) -> str:
    words = value.replace("-", " ").replace("_", " ").split()
    acronyms = {"mpi": "MPI", "openmpi": "OpenMPI", "mpich": "MPICH"}
    return " ".join(acronyms.get(word.lower(), word.title()) for word in words)


def _profile_id(target: ImageTarget) -> str:
    return f"{target.flavor}-{target.dist}".replace(":", "-")


def _profile_name(target: ImageTarget) -> str:
    if target.flavor == "spack":
        return f"Feel++ Spack {_title(target.dist)}"
    suffix = f" ({target.version})" if target.version else ""
    return f"Feel++ {_title(target.flavor)} {_title(target.dist)}{suffix}"


def _base_cmake_preset(target: ImageTarget) -> str:
    if target.metadata.get("cmake_devcontainer_preset"):
        return target.metadata["cmake_devcontainer_preset"]
    if target.image_backend == "spack":
        return target.metadata.get("cmake_full_preset", "release-clang-spack")
    return "default"


def _profile_from_target(target: ImageTarget) -> DevcontainerProfile:
    profile_id = _profile_id(target)
    return DevcontainerProfile(
        id=profile_id,
        target=target.target,
        name=_profile_name(target),
        image=oci_image_ref("feelpp-env", target.oci_dist),
        cmake_preset=f"container-{profile_id}",
        inherits_preset=_base_cmake_preset(target),
        distribution_env=f"-{target.oci_dist}",
        build_dir=f"${{sourceDir}}/build/container/{profile_id}",
        install_dir=f"${{sourceDir}}/install/container/{profile_id}",
    )


def default_targets(repo_root: str | Path | None = None) -> list[str]:
    root = _repo_root(repo_root)
    profile = load_image_profile_config(root, profile=DEFAULT_IMAGE_PROFILE)
    groups = profile.get("groups", {})
    if isinstance(groups, dict) and DEFAULT_DEVCONTAINER_GROUP in groups:
        targets = groups[DEFAULT_DEVCONTAINER_GROUP]
        if isinstance(targets, list):
            return [str(target) for target in targets]
    defaults = profile.get("defaults", {})
    if isinstance(defaults, dict) and isinstance(defaults.get("targets"), list):
        return [str(target) for target in defaults["targets"]]
    return [DEFAULT_DEVCONTAINER_TARGET]


def load_profiles(
    repo_root: str | Path | None = None,
    *,
    targets: Iterable[str] | None = None,
) -> list[DevcontainerProfile]:
    root = _repo_root(repo_root)
    selected_targets = list(targets) if targets is not None else default_targets(root)
    profiles = [
        _profile_from_target(
            get_image_target(target, repo_root=root, profile=DEFAULT_IMAGE_PROFILE)
        )
        for target in selected_targets
    ]
    return sorted(profiles, key=lambda profile: profile.id)


def default_profile(
    profiles: Iterable[DevcontainerProfile],
    *,
    default_target: str = DEFAULT_DEVCONTAINER_TARGET,
) -> DevcontainerProfile:
    rows = list(profiles)
    for profile in rows:
        if profile.target == default_target:
            return profile
    if not rows:
        raise ValueError("No devcontainer profiles selected")
    return rows[0]


def render_devcontainer(profile: DevcontainerProfile, *, workspace_name: str) -> dict[str, object]:
    return {
        "name": profile.name,
        "image": profile.image,
        "workspaceFolder": f"/workspaces/{workspace_name}",
        "remoteUser": profile.remote_user,
        "updateRemoteUserUID": True,
        "init": True,
        "containerEnv": {
            "FEELPP_CMAKE_PRESET": profile.cmake_preset,
            "DISTRIBUTION": profile.distribution_env,
        },
        "remoteEnv": {
            "CMAKE_BUILD_PARALLEL_LEVEL": str(profile.build_jobs),
        },
        "postCreateCommand": (
            "${containerWorkspaceFolder}/.devcontainer/feelppconfig.sh "
            "${containerWorkspaceFolder}"
        ),
        "customizations": {
            "vscode": {
                "settings": {
                    "cmake.configurePreset": profile.cmake_preset,
                    "cmake.buildPreset": profile.cmake_preset,
                },
                "extensions": VSCODE_EXTENSIONS,
            }
        },
    }


def render_cmake_user_presets(profiles: Iterable[DevcontainerProfile]) -> dict[str, object]:
    rows = list(profiles)
    return {
        "version": 6,
        "configurePresets": [
            {
                "name": profile.cmake_preset,
                "displayName": f"{profile.name} container",
                "description": (
                    f"{profile.name} dev container configuration using a dedicated "
                    "container build directory."
                ),
                "inherits": profile.inherits_preset,
                "binaryDir": profile.build_dir,
                "cacheVariables": {
                    "CMAKE_INSTALL_PREFIX": profile.install_dir,
                },
            }
            for profile in rows
        ],
        "buildPresets": [
            {
                "name": profile.cmake_preset,
                "configurePreset": profile.cmake_preset,
                "jobs": profile.build_jobs,
            }
            for profile in rows
        ],
    }


def render_files(
    repo_root: str | Path | None = None,
    *,
    targets: Iterable[str] | None = None,
    default_target: str = DEFAULT_DEVCONTAINER_TARGET,
) -> dict[Path, str]:
    root = _repo_root(repo_root)
    profiles = load_profiles(root, targets=targets)
    workspace_name = root.name
    files: dict[Path, str] = {}
    for profile in profiles:
        files[profile.relative_path] = _json(
            render_devcontainer(profile, workspace_name=workspace_name)
        )
    files[ROOT_DEVCONTAINER_PATH] = _json(
        render_devcontainer(
            default_profile(profiles, default_target=default_target),
            workspace_name=workspace_name,
        )
    )
    files[CMAKE_USER_PRESETS_PATH] = _json(render_cmake_user_presets(profiles))
    return files


@dataclass(frozen=True)
class WriteResult:
    written: list[Path]
    unchanged: list[Path]
    drifted: list[Path]


def write_generated_files(
    repo_root: str | Path | None = None,
    *,
    targets: Iterable[str] | None = None,
    default_target: str = DEFAULT_DEVCONTAINER_TARGET,
    check: bool = False,
) -> WriteResult:
    root = _repo_root(repo_root)
    rendered = render_files(root, targets=targets, default_target=default_target)
    written: list[Path] = []
    unchanged: list[Path] = []
    drifted: list[Path] = []
    for relative_path, contents in rendered.items():
        path = root / relative_path
        current = path.read_text(encoding="utf-8") if path.is_file() else None
        if current == contents:
            unchanged.append(relative_path)
            continue
        if check:
            drifted.append(relative_path)
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents, encoding="utf-8")
        written.append(relative_path)
    return WriteResult(written=written, unchanged=unchanged, drifted=drifted)


def validate_generated_files(
    repo_root: str | Path | None = None,
    *,
    targets: Iterable[str] | None = None,
    default_target: str = DEFAULT_DEVCONTAINER_TARGET,
) -> WriteResult:
    return write_generated_files(
        repo_root,
        targets=targets,
        default_target=default_target,
        check=True,
    )
