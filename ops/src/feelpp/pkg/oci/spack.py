from __future__ import annotations

import json
from pathlib import Path
import shutil

from .catalog import ImageTarget, list_image_targets
from .common import default_bake_target_name, oci_image_ref
from ..core.context import WorkspaceContext


DEFAULT_SPACK_REF = "v1.0.0"
SHARED_SPACK_DIR = Path("packaging") / "spack"


def _environment_manifest(repo_root: Path, environment_name: str) -> Path:
    manifest = repo_root / "packaging" / "spack" / "environments" / environment_name / "spack.yaml"
    if not manifest.is_file():
        raise ValueError(
            f"Unsupported Spack environment '{environment_name}': {manifest} does not exist"
        )
    return manifest


def _ignore_generated_spack_state(_root: str, names: list[str]) -> set[str]:
    ignored = {".spack-env", "__pycache__"}
    return {name for name in names if name in ignored}


def render_spack_dockerfile(
    *,
    base_image: str,
    spack_ref: str,
    environment_name: str,
) -> str:
    environment_dir = f"/opt/feelpp/packaging/spack/environments/{environment_name}"
    return f"""# syntax=docker/dockerfile:1
ARG BASE_IMAGE={base_image}
FROM ${{BASE_IMAGE}}

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \\
    && apt-get install -y --no-install-recommends \\
       bzip2 \\
       ca-certificates \\
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
RUN git clone --branch "${{SPACK_REF}}" --depth=1 https://github.com/spack/spack.git /opt/spack

ENV SPACK_ROOT=/opt/spack \\
    SPACK_USER_CONFIG_PATH=/opt/spack-user \\
    SPACK_USER_CACHE_PATH=/opt/spack-user-cache

COPY packaging/spack /opt/feelpp/packaging/spack

RUN mkdir -p "$SPACK_USER_CONFIG_PATH" "$SPACK_USER_CACHE_PATH" \\
    && . "$SPACK_ROOT/share/spack/setup-env.sh" \\
    && spack -e {environment_dir} concretize \\
    && spack -e {environment_dir} install \\
    && spack clean --all

RUN cat >/etc/profile.d/feelpp-spack.sh <<'EOF'
export SPACK_ROOT=/opt/spack
export SPACK_USER_CONFIG_PATH=/opt/spack-user
export SPACK_USER_CACHE_PATH=/opt/spack-user-cache
. "$SPACK_ROOT/share/spack/setup-env.sh"
spack env activate {environment_dir}
EOF

CMD ["/bin/bash"]
"""


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
) -> dict[str, object]:
    resolved_environment = environment_name or target.spack_environment or f"cpu/{target.dist}"
    environment_manifest = _environment_manifest(workspace.repo_root, resolved_environment)
    resolved_base_image = base_image or target.base_image
    resolved_bake_target = bake_target or default_bake_target_name(target.target)
    resolved_image_tag = image_tag or oci_image_ref(
        "feelpp-env",
        target.oci_dist,
        registry=registry,
        namespace=namespace,
    )

    context_dir = workspace.job_root / "images" / resolved_bake_target
    packaging_dst = context_dir / SHARED_SPACK_DIR
    shutil.rmtree(context_dir, ignore_errors=True)
    packaging_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        workspace.repo_root / SHARED_SPACK_DIR,
        packaging_dst,
        ignore=_ignore_generated_spack_state,
    )

    dockerfile_path = context_dir / "Dockerfile"
    dockerfile_path.write_text(
        render_spack_dockerfile(
            base_image=resolved_base_image,
            spack_ref=spack_ref,
            environment_name=resolved_environment,
        ),
        encoding="utf-8",
    )

    bake_payload = {
        "group": {
            "default": {
                "targets": [resolved_bake_target],
            }
        },
        "target": {
            resolved_bake_target: {
                "context": ".",
                "dockerfile": "Dockerfile",
                "tags": [resolved_image_tag],
                "args": {
                    "BASE_IMAGE": resolved_base_image,
                    "SPACK_REF": spack_ref,
                },
            }
        }
    }
    bake_file = context_dir / "docker-bake.json"
    bake_file.write_text(json.dumps(bake_payload, indent=2) + "\n", encoding="utf-8")
    recommended_groups = ["default"]

    return {
        "repo_root": str(workspace.repo_root),
        "job_root": str(workspace.job_root),
        "packaging_target": target.target,
        "image_backend": target.image_backend,
        "image_strategy": target.image_strategy,
        "bake_target": resolved_bake_target,
        "oci_dist": target.oci_dist,
        "image_tag": resolved_image_tag,
        "base_image": resolved_base_image,
        "spack_ref": spack_ref,
        "environment": resolved_environment,
        "environment_manifest": str(environment_manifest),
        "context_dir": str(context_dir),
        "dockerfile": str(dockerfile_path),
        "bake_file": str(bake_file),
        "default_group": "default",
        "available_groups": ["default"],
        "recommended_groups": recommended_groups,
        "docker_bake_command": f"docker buildx bake -f {bake_file} {' '.join(recommended_groups)}",
    }
