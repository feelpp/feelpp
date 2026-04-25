from __future__ import annotations

import json
from pathlib import Path
import shutil

import yaml

from .bake_permissions import fs_read_allow_flags, required_fs_read_paths_from_bake_payload
from .catalog import ImageTarget, list_image_targets
from .cmake_presets import resolve_cmake_preset
from .common import branch_tag_suffix, default_bake_target_name, oci_image_ref
from .platforms import normalize_platform_overrides
from ..core.context import WorkspaceContext


DEFAULT_SPACK_REF = "v1.0.0"
DEFAULT_SPACK_BUILD_JOBS = 16
DEFAULT_SPACK_CONCURRENT_PACKAGES = 0
DEFAULT_SPACK_FAIL_FAST = True
DEFAULT_SPACK_SHOW_LOG_ON_ERROR = True
SHARED_SPACK_DIR = Path("packaging") / "spack"
DOCKER_TEMPLATE_DIR = Path("packaging") / "docker" / "templates"


def _environment_manifest(repo_root: Path, environment_name: str) -> Path:
    manifest = repo_root / "packaging" / "spack" / "environments" / environment_name / "spack.yaml"
    if not manifest.is_file():
        raise ValueError(
            f"Unsupported Spack environment '{environment_name}': {manifest} does not exist"
        )
    return manifest


def _ignore_generated_spack_state(_root: str, names: list[str]) -> set[str]:
    ignored = {".spack-env", "__pycache__", "spack.lock"}
    return {name for name in names if name in ignored}


def _read_spack_manifest(path: Path) -> dict[str, object]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid Spack manifest: {path}")
    return payload


def _resolved_parallelism(
    manifest_path: Path,
    *,
    build_jobs: int | None = None,
    concurrent_packages: int | None = None,
) -> tuple[int, int]:
    payload = _read_spack_manifest(manifest_path)
    spack_payload = payload.get("spack", {}) or {}
    if not isinstance(spack_payload, dict):
        raise ValueError(f"Invalid spack section in {manifest_path}")
    config_payload = spack_payload.get("config", {}) or {}
    if not isinstance(config_payload, dict):
        raise ValueError(f"Invalid spack.config section in {manifest_path}")

    resolved_build_jobs = (
        build_jobs
        if build_jobs is not None
        else int(config_payload.get("build_jobs", DEFAULT_SPACK_BUILD_JOBS))
    )
    resolved_concurrent_packages = (
        concurrent_packages
        if concurrent_packages is not None
        else int(config_payload.get("concurrent_packages", DEFAULT_SPACK_CONCURRENT_PACKAGES))
    )

    if resolved_build_jobs < 1:
        raise ValueError("Spack build_jobs must be >= 1")
    if resolved_concurrent_packages < 0:
        raise ValueError("Spack concurrent_packages must be >= 0")

    return resolved_build_jobs, resolved_concurrent_packages


def render_spack_dockerfile(
    *,
    base_image: str,
    spack_ref: str,
    environment_name: str,
    build_jobs: int,
    concurrent_packages: int,
    fail_fast: bool,
    show_log_on_error: bool,
) -> str:
    environment_dir = f"/opt/feelpp/packaging/spack/environments/{environment_name}"
    install_flag_parts: list[str] = []
    if fail_fast:
        install_flag_parts.append("--fail-fast")
    if show_log_on_error:
        install_flag_parts.append("--show-log-on-error")
    install_clause = " ".join(["install", *install_flag_parts])
    return f"""# syntax=docker/dockerfile:1
ARG BASE_IMAGE={base_image}
FROM ${{BASE_IMAGE}}

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \\
    && apt-get install -y --no-install-recommends \\
       bzip2 \\
       ca-certificates \\
       clang \\
       file \\
       g++ \\
       gcc \\
       gfortran \\
       git \\
       gzip \\
       patch \\
       python3 \\
       tar \\
       unzip \\
       xz-utils \\
       zstd \\
    && rm -rf /var/lib/apt/lists/*

ARG SPACK_REF={spack_ref}
ARG SPACK_BUILD_JOBS={build_jobs}
ARG SPACK_CONCURRENT_PACKAGES={concurrent_packages}
RUN git clone --branch "${{SPACK_REF}}" --depth=1 https://github.com/spack/spack.git /opt/spack

ENV SPACK_ROOT=/opt/spack \\
    SPACK_USER_CONFIG_PATH=/opt/spack-user \\
    SPACK_USER_CACHE_PATH=/opt/spack-user-cache

COPY packaging/spack /opt/feelpp/packaging/spack

RUN mkdir -p "$SPACK_USER_CONFIG_PATH" "$SPACK_USER_CACHE_PATH" \\
    && . "$SPACK_ROOT/share/spack/setup-env.sh" \\
    && spack -e {environment_dir} concretize -f \\
    && if [ "${{SPACK_CONCURRENT_PACKAGES}}" -gt 0 ]; then \\
         spack -e {environment_dir} {install_clause} -j "${{SPACK_BUILD_JOBS}}" -p "${{SPACK_CONCURRENT_PACKAGES}}"; \\
       else \\
         spack -e {environment_dir} {install_clause} -j "${{SPACK_BUILD_JOBS}}"; \\
       fi \\
    && spack clean --all

RUN cat >/etc/profile.d/feelpp-spack.sh <<'EOF'
export SPACK_ROOT=/opt/spack
export SPACK_USER_CONFIG_PATH=/opt/spack-user
export SPACK_USER_CACHE_PATH=/opt/spack-user-cache
. "$SPACK_ROOT/share/spack/setup-env.sh"
spack env activate {environment_dir}
EOF

ENV BASH_ENV=/etc/profile.d/feelpp-spack.sh
SHELL ["/bin/bash", "-lc"]

CMD ["/bin/bash"]
"""


def _normalize_requested_component(component: str | None) -> str | None:
    if component in {None, ""}:
        return None
    normalized = str(component).strip().lower()
    if normalized in {"env"}:
        return "env"
    if normalized == "full":
        return "full"
    raise ValueError("Spack image generation only supports component values env or full")


def _full_builder_tag(target: ImageTarget, *, branch: str) -> str:
    return f"{target.oci_dist}{branch_tag_suffix(branch)}-full-dev"


def _full_runtime_tag(target: ImageTarget, *, branch: str) -> str:
    return f"{target.oci_dist}{branch_tag_suffix(branch)}-full"


def describe_spack_target(repo_root: Path, target: ImageTarget) -> dict[str, object]:
    environment_name = target.spack_environment or f"cpu/{target.dist}"
    return {
        **target.as_dict(),
        "environment": environment_name,
        "environment_manifest": str(
            repo_root / "packaging" / "spack" / "environments" / environment_name / "spack.yaml"
        ),
        "supported": _environment_manifest(repo_root, environment_name).is_file(),
    }


def list_spack_targets(repo_root: Path) -> list[ImageTarget]:
    return list_image_targets(repo_root, backend="spack")


def generate_spack_bake(
    *,
    workspace: WorkspaceContext,
    target: ImageTarget,
    environment_name: str | None = None,
    base_image: str | None = None,
    spack_ref: str = DEFAULT_SPACK_REF,
    bake_target: str | None = None,
    image_tag: str | None = None,
    registry: str | None = None,
    namespace: str | None = None,
    component: str | None = None,
    from_image: str | None = None,
    cc: str = "clang",
    cxx: str = "clang++",
    cmake_flags: str = "",
    spack_build_jobs: int | None = None,
    spack_concurrent_packages: int | None = None,
    spack_fail_fast: bool | None = None,
    spack_show_log_on_error: bool | None = None,
    platform_overrides: list[str] | None = None,
) -> dict[str, object]:
    resolved_environment = environment_name or target.spack_environment or f"cpu/{target.dist}"
    environment_manifest = _environment_manifest(workspace.repo_root, resolved_environment)
    resolved_base_image = base_image or target.base_image
    resolved_bake_target = bake_target or default_bake_target_name(target.target)
    requested_component = _normalize_requested_component(component)
    resolved_build_jobs, resolved_concurrent_packages = _resolved_parallelism(
        environment_manifest,
        build_jobs=spack_build_jobs,
        concurrent_packages=spack_concurrent_packages,
    )
    resolved_fail_fast = (
        DEFAULT_SPACK_FAIL_FAST if spack_fail_fast is None else bool(spack_fail_fast)
    )
    resolved_show_log_on_error = (
        DEFAULT_SPACK_SHOW_LOG_ON_ERROR
        if spack_show_log_on_error is None
        else bool(spack_show_log_on_error)
    )
    resolved_platforms = normalize_platform_overrides(platform_overrides)
    env_image_ref = image_tag or oci_image_ref(
        "feelpp-env",
        target.oci_dist,
        registry=registry,
        namespace=namespace,
    )

    context_dir = workspace.job_root / "images" / resolved_bake_target
    shutil.rmtree(context_dir, ignore_errors=True)
    context_dir.mkdir(parents=True, exist_ok=True)

    dockerfile_path = context_dir / "Dockerfile"
    available_groups = ["default"]
    recommended_groups = ["default"]
    image_refs: dict[str, str] = {"feelpp-env": env_image_ref}

    if requested_component in {None, "env"}:
        packaging_dst = context_dir / SHARED_SPACK_DIR
        packaging_dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(
            workspace.repo_root / SHARED_SPACK_DIR,
            packaging_dst,
            ignore=_ignore_generated_spack_state,
        )

        dockerfile_path.write_text(
            render_spack_dockerfile(
                base_image=resolved_base_image,
                spack_ref=spack_ref,
                environment_name=resolved_environment,
                build_jobs=resolved_build_jobs,
                concurrent_packages=resolved_concurrent_packages,
                fail_fast=resolved_fail_fast,
                show_log_on_error=resolved_show_log_on_error,
            ),
            encoding="utf-8",
        )

        target_payload = {
            "context": ".",
            "dockerfile": "Dockerfile",
            "tags": [env_image_ref],
            "args": {
                "BASE_IMAGE": resolved_base_image,
                "SPACK_REF": spack_ref,
                "SPACK_BUILD_JOBS": str(resolved_build_jobs),
                "SPACK_CONCURRENT_PACKAGES": str(resolved_concurrent_packages),
            },
        }
        if resolved_platforms is not None:
            target_payload["platforms"] = list(resolved_platforms)

        bake_payload = {
            "group": {
                "default": {
                    "targets": [resolved_bake_target],
                }
            },
            "target": {
                resolved_bake_target: target_payload,
            }
        }
    else:
        feelpp_dir = context_dir / "feelpp"
        feelpp_dir.mkdir(parents=True, exist_ok=True)
        full_template = workspace.repo_root / DOCKER_TEMPLATE_DIR / "feelpp.Dockerfile.full"
        dockerfile_path = feelpp_dir / "Dockerfile.full"
        dockerfile_path.write_text(full_template.read_text(encoding="utf-8"), encoding="utf-8")

        full_builder_ref = oci_image_ref(
            "feelpp",
            _full_builder_tag(target, branch=workspace.branch),
            registry=registry,
            namespace=namespace,
        )
        full_runtime_ref = oci_image_ref(
            "feelpp",
            _full_runtime_tag(target, branch=workspace.branch),
            registry=registry,
            namespace=namespace,
        )
        image_refs["feelpp-full-dev"] = full_builder_ref
        image_refs["feelpp-full"] = full_runtime_ref
        recommended_groups = ["full-all"]
        available_groups = ["full-dev", "full-runtime", "full-all"]

        full_builder_target = {
            "context": "feelpp",
            "dockerfile": "Dockerfile.full",
            "target": "builder",
            "contexts": {"feelpp_source": str(workspace.repo_root)},
            "args": {
                "FROM_IMAGE": from_image or env_image_ref,
                "DESCRIPTION": "Feel++ Full Stack (dev)",
                "BRANCH": workspace.branch,
                "CMAKE_PRESET": resolve_cmake_preset(target, component="full"),
                "CXX": cxx,
                "CC": cc,
                "CMAKE_FLAGS": cmake_flags,
            },
            "tags": [full_builder_ref],
        }
        full_runtime_target = {
            "context": "feelpp",
            "dockerfile": "Dockerfile.full",
            "target": "runtime",
            "contexts": {
                "feelpp_source": str(workspace.repo_root),
                "build_output": "target:feelpp-full",
            },
            "args": {
                "FROM_IMAGE": from_image or env_image_ref,
                "DESCRIPTION": "Feel++ Full Stack",
                "BRANCH": workspace.branch,
                "CMAKE_PRESET": resolve_cmake_preset(target, component="full"),
                "CXX": cxx,
                "CC": cc,
                "CMAKE_FLAGS": cmake_flags,
            },
            "tags": [full_runtime_ref],
        }
        if resolved_platforms is not None:
            full_builder_target["platforms"] = list(resolved_platforms)
            full_runtime_target["platforms"] = list(resolved_platforms)

        bake_payload = {
            "group": {
                "full-dev": {"targets": ["feelpp-full"]},
                "full-runtime": {"targets": ["feelpp-full-runtime"]},
                "full-all": {"targets": ["feelpp-full", "feelpp-full-runtime"]},
                "default": {"targets": ["feelpp-full", "feelpp-full-runtime"]},
            },
            "target": {
                "feelpp-full": full_builder_target,
                "feelpp-full-runtime": full_runtime_target,
            },
        }
    bake_file = context_dir / "docker-bake.json"
    bake_file.write_text(json.dumps(bake_payload, indent=2) + "\n", encoding="utf-8")
    fs_read_paths = required_fs_read_paths_from_bake_payload(
        bake_payload,
        base_dir=context_dir,
        groups=recommended_groups,
    )
    docker_bake_command = ["docker", "buildx", "bake", *fs_read_allow_flags(fs_read_paths), "-f", str(bake_file)]
    docker_bake_command.extend(recommended_groups)

    return {
        "repo_root": str(workspace.repo_root),
        "job_root": str(workspace.job_root),
        "packaging_target": target.target,
        "image_backend": target.image_backend,
        "image_strategy": target.image_strategy,
        "bake_target": resolved_bake_target,
        "oci_dist": target.oci_dist,
        "image_tag": env_image_ref,
        "base_image": resolved_base_image,
        "spack_ref": spack_ref,
        "spack_build_jobs": resolved_build_jobs,
        "spack_concurrent_packages": resolved_concurrent_packages,
        "spack_fail_fast": resolved_fail_fast,
        "spack_show_log_on_error": resolved_show_log_on_error,
        "selected_component": requested_component or "env",
        "from_image": from_image or "",
        "platforms": list(resolved_platforms or []),
        "environment": resolved_environment,
        "environment_manifest": str(environment_manifest),
        "context_dir": str(context_dir),
        "dockerfile": str(dockerfile_path),
        "bake_file": str(bake_file),
        "default_group": "default",
        "available_groups": available_groups,
        "recommended_groups": recommended_groups,
        "fs_read_paths": fs_read_paths,
        "image_refs": image_refs,
        "docker_bake_command": " ".join(docker_bake_command),
    }
