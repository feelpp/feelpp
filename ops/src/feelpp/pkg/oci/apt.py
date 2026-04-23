from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import shutil

from jinja2 import Environment, FileSystemLoader, StrictUndefined

from .catalog import ImageTarget, list_image_targets
from .common import (
    branch_tag_suffix,
    default_bake_target_name,
    oci_image_ref,
)
from .docker_metadata import (
    DockerMetadata,
    load_docker_metadata,
    packages_for_variant,
    resolve_distribution_version,
)
from ..core.context import WorkspaceContext
from ..graph import load_manifest


DEFAULT_VARIANT = "feelpp-env"
DEFAULT_CC = "clang"
DEFAULT_CXX = "clang++"
SHARED_ASSET_NAMES = (
    "WELCOME",
    "bashrc.feelpp",
    "ctest-to-junit.xsl",
    "feelpp.conf.sh",
    "feelpp.env.sh",
    "start-user.sh",
    "start.sh",
)


@dataclass(frozen=True)
class ComponentTemplateSpec:
    component_name: str
    bake_target: str
    runtime_target: str
    image_name: str
    template_name: str
    description: str


COMPONENT_TEMPLATE_SPECS = {
    "feelpp": ComponentTemplateSpec(
        component_name="feelpp",
        bake_target="feelpp",
        runtime_target="feelpp-runtime",
        image_name="feelpp",
        template_name="feelpp.Dockerfile.multistage",
        description="Feel++ core library",
    ),
    "feelpp-toolboxes": ComponentTemplateSpec(
        component_name="feelpp-toolboxes",
        bake_target="toolboxes",
        runtime_target="toolboxes-runtime",
        image_name="feelpp-toolboxes",
        template_name="feelpp-toolboxes.Dockerfile.multistage",
        description="Feel++ Toolboxes",
    ),
    "feelpp-mor": ComponentTemplateSpec(
        component_name="feelpp-mor",
        bake_target="mor",
        runtime_target="mor-runtime",
        image_name="feelpp-mor",
        template_name="feelpp-mor.Dockerfile.multistage",
        description="Feel++ Model Order Reduction",
    ),
}


def _join_packages(packages: list[str], indent: int = 8) -> str:
    spaces = " " * indent
    return "\n".join(f"{spaces}{package} \\" for package in packages)


def _template_environment(metadata: DockerMetadata) -> Environment:
    env = Environment(
        loader=FileSystemLoader(str(metadata.root / "templates")),
        undefined=StrictUndefined,
        trim_blocks=True,
        lstrip_blocks=True,
    )
    env.filters["join_packages"] = _join_packages
    return env


def _component_template_path(metadata: DockerMetadata, template_name: str) -> Path:
    return metadata.root / "templates" / template_name


def _copy_shared_assets(metadata: DockerMetadata, destination: Path) -> None:
    assets_dir = metadata.root / "assets"
    for asset_name in SHARED_ASSET_NAMES:
        source = assets_dir / asset_name
        if not source.is_file():
            raise ValueError(f"Missing Docker shared asset: {source}")
        shutil.copyfile(source, destination / asset_name)


def _copy_template(source: Path, destination: Path) -> None:
    destination.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")


def _component_set(component_set_name: str | None, repo_root: Path) -> list[str]:
    if component_set_name not in {None, "", "default"}:
        raise ValueError(f"Unsupported component set: {component_set_name}")
    manifest = load_manifest(repo_root / "packaging" / "manifest" / "components.toml")
    return list(manifest.default_components)


def _builder_tag(target: ImageTarget, *, branch: str) -> str:
    return f"{target.oci_dist}{branch_tag_suffix(branch)}-dev"


def _runtime_tag(target: ImageTarget, *, branch: str) -> str:
    return f"{target.oci_dist}{branch_tag_suffix(branch)}"


def _full_builder_tag(target: ImageTarget, *, branch: str) -> str:
    return f"{target.oci_dist}{branch_tag_suffix(branch)}-full-dev"


def _full_runtime_tag(target: ImageTarget, *, branch: str) -> str:
    return f"{target.oci_dist}{branch_tag_suffix(branch)}-full"


def _environment_image_ref(
    target: ImageTarget,
    *,
    registry: str | None = None,
    namespace: str | None = None,
) -> str:
    return oci_image_ref("feelpp-env", target.oci_dist, registry=registry, namespace=namespace)


def _component_image_ref(
    image_name: str,
    tag: str,
    *,
    registry: str | None = None,
    namespace: str | None = None,
) -> str:
    return oci_image_ref(image_name, tag, registry=registry, namespace=namespace)


def describe_apt_target(
    repo_root: Path,
    target: ImageTarget,
    *,
    variant_name: str = DEFAULT_VARIANT,
) -> dict[str, object]:
    metadata = load_docker_metadata(repo_root)
    distribution = resolve_distribution_version(metadata, family=target.flavor, name=target.dist)
    packages = packages_for_variant(metadata, variant_name=variant_name)
    return {
        **target.as_dict(),
        "variant": variant_name,
        "platforms": list(distribution.platforms),
        "template": distribution.template,
        "packages_count": len(packages),
        "supported": True,
    }


def list_apt_targets(repo_root: Path) -> list[ImageTarget]:
    return list_image_targets(repo_root, backend="apt")


def _render_environment_dockerfile(
    *,
    metadata: DockerMetadata,
    target: ImageTarget,
    variant_name: str,
    base_image: str | None = None,
) -> str:
    distribution = resolve_distribution_version(metadata, family=target.flavor, name=target.dist)
    variant = metadata.variants[variant_name]
    template = _template_environment(metadata).get_template(distribution.template)
    packages = packages_for_variant(metadata, variant_name=variant_name)
    openmpi_packages = list(distribution.openmpi_packages) if variant.enable_openmpi else []
    return template.render(
        base_image=base_image or distribution.base_image or target.base_image,
        dist_name=distribution.family,
        dist_version=distribution.version,
        dist_codename=distribution.name,
        variant_name=variant.description or variant.name,
        packages=packages,
        openmpi_packages=openmpi_packages,
    )


def _component_targets_for_row(repo_root: Path, target: ImageTarget) -> list[ComponentTemplateSpec]:
    component_names = _component_set(target.component_set, repo_root)
    component_specs: list[ComponentTemplateSpec] = []
    for component_name in component_names:
        try:
            component_specs.append(COMPONENT_TEMPLATE_SPECS[component_name])
        except KeyError as exc:
            raise ValueError(f"No Docker component template mapping for {component_name}") from exc
    return component_specs


def generate_apt_bake(
    *,
    workspace: WorkspaceContext,
    target: ImageTarget,
    variant_name: str = DEFAULT_VARIANT,
    cc: str = DEFAULT_CC,
    cxx: str = DEFAULT_CXX,
    cmake_flags: str = "",
    base_image: str | None = None,
    registry: str | None = None,
    namespace: str | None = None,
    bake_target: str | None = None,
) -> dict[str, object]:
    metadata = load_docker_metadata(workspace.repo_root)
    distribution = resolve_distribution_version(metadata, family=target.flavor, name=target.dist)
    component_specs = _component_targets_for_row(workspace.repo_root, target)

    resolved_bake_target = bake_target or default_bake_target_name(target.target)
    context_dir = workspace.job_root / "images" / resolved_bake_target
    env_dir = context_dir / "feelpp-env"
    feelpp_dir = context_dir / "feelpp"
    toolboxes_dir = context_dir / "feelpp-toolboxes"
    mor_dir = context_dir / "feelpp-mor"

    shutil.rmtree(context_dir, ignore_errors=True)
    for directory in (env_dir, feelpp_dir, toolboxes_dir, mor_dir):
        directory.mkdir(parents=True, exist_ok=True)

    env_dockerfile = env_dir / "Dockerfile"
    env_dockerfile.write_text(
        _render_environment_dockerfile(
            metadata=metadata,
            target=target,
            variant_name=variant_name,
            base_image=base_image or target.base_image,
        ),
        encoding="utf-8",
    )
    _copy_shared_assets(metadata, env_dir)

    _copy_template(
        _component_template_path(metadata, "feelpp.Dockerfile.multistage"),
        feelpp_dir / "Dockerfile.multistage",
    )
    _copy_template(
        _component_template_path(metadata, "feelpp.Dockerfile.full"),
        feelpp_dir / "Dockerfile.full",
    )
    _copy_template(
        _component_template_path(metadata, "feelpp-toolboxes.Dockerfile.multistage"),
        toolboxes_dir / "Dockerfile.multistage",
    )
    _copy_template(
        _component_template_path(metadata, "feelpp-mor.Dockerfile.multistage"),
        mor_dir / "Dockerfile.multistage",
    )

    targets: dict[str, dict[str, object]] = {
        "feelpp-env": {
            "context": "feelpp-env",
            "dockerfile": "Dockerfile",
            "tags": [_environment_image_ref(target, registry=registry, namespace=namespace)],
            "platforms": list(distribution.platforms),
        }
    }

    previous_runtime_target = "feelpp-env"
    previous_context_alias = "feelpp_env_image"
    for component_spec in component_specs:
        context_dir_name = {
            "feelpp": "feelpp",
            "feelpp-toolboxes": "feelpp-toolboxes",
            "feelpp-mor": "feelpp-mor",
        }[component_spec.component_name]
        builder_tag = _builder_tag(target, branch=workspace.branch)
        runtime_tag = _runtime_tag(target, branch=workspace.branch)
        if component_spec.component_name == "feelpp":
            builder_image_name = "feelpp"
            runtime_image_name = "feelpp"
        elif component_spec.component_name == "feelpp-toolboxes":
            builder_image_name = "feelpp-toolboxes"
            runtime_image_name = "feelpp-toolboxes"
        else:
            builder_image_name = "feelpp-mor"
            runtime_image_name = "feelpp-mor"

        contexts = {previous_context_alias: f"target:{previous_runtime_target}"}
        targets[component_spec.bake_target] = {
            "context": context_dir_name,
            "dockerfile": "Dockerfile.multistage",
            "target": "builder",
            "contexts": contexts,
            "args": {
                "FROM_IMAGE": previous_context_alias,
                "DESCRIPTION": f"{component_spec.description} (dev)",
                "BRANCH": workspace.branch,
                "CXX": cxx,
                "CC": cc,
                "CMAKE_FLAGS": cmake_flags,
            },
            "tags": [
                _component_image_ref(
                    builder_image_name,
                    builder_tag,
                    registry=registry,
                    namespace=namespace,
                )
            ],
            "platforms": list(distribution.platforms),
        }
        targets[component_spec.runtime_target] = {
            "context": context_dir_name,
            "dockerfile": "Dockerfile.multistage",
            "target": "runtime",
            "contexts": contexts,
            "args": {
                "FROM_IMAGE": previous_context_alias,
                "DESCRIPTION": component_spec.description,
                "BRANCH": workspace.branch,
                "CXX": cxx,
                "CC": cc,
                "CMAKE_FLAGS": cmake_flags,
            },
            "tags": [
                _component_image_ref(
                    runtime_image_name,
                    runtime_tag,
                    registry=registry,
                    namespace=namespace,
                )
            ],
            "platforms": list(distribution.platforms),
        }
        previous_runtime_target = component_spec.runtime_target
        previous_context_alias = f"{component_spec.bake_target}_runtime_image"

    targets["feelpp-full"] = {
        "context": "feelpp",
        "dockerfile": "Dockerfile.full",
        "target": "builder",
        "contexts": {"feelpp_env_image": "target:feelpp-env"},
        "args": {
            "FROM_IMAGE": "feelpp_env_image",
            "DESCRIPTION": "Feel++ Full Stack (dev)",
            "BRANCH": workspace.branch,
            "CXX": cxx,
            "CC": cc,
            "CMAKE_FLAGS": cmake_flags,
        },
        "tags": [
            _component_image_ref(
                "feelpp",
                _full_builder_tag(target, branch=workspace.branch),
                registry=registry,
                namespace=namespace,
            )
        ],
        "platforms": list(distribution.platforms),
    }
    targets["feelpp-full-runtime"] = {
        "context": "feelpp",
        "dockerfile": "Dockerfile.full",
        "target": "runtime",
        "contexts": {"feelpp_env_image": "target:feelpp-env"},
        "args": {
            "FROM_IMAGE": "feelpp_env_image",
            "DESCRIPTION": "Feel++ Full Stack",
            "BRANCH": workspace.branch,
            "CXX": cxx,
            "CC": cc,
            "CMAKE_FLAGS": cmake_flags,
        },
        "tags": [
            _component_image_ref(
                "feelpp",
                _full_runtime_tag(target, branch=workspace.branch),
                registry=registry,
                namespace=namespace,
            )
        ],
        "platforms": list(distribution.platforms),
    }

    bake_payload = {
        "group": {
            "default": {"targets": ["feelpp"]},
            "all-dev": {"targets": [spec.bake_target for spec in component_specs]},
            "all": {"targets": [spec.runtime_target for spec in component_specs]},
            "full": {"targets": ["feelpp-full"]},
            "full-all": {"targets": ["feelpp-full", "feelpp-full-runtime"]},
        },
        "target": targets,
    }
    bake_file = context_dir / "docker-bake.json"
    bake_file.write_text(json.dumps(bake_payload, indent=2) + "\n", encoding="utf-8")

    return {
        "repo_root": str(workspace.repo_root),
        "job_root": str(workspace.job_root),
        "packaging_target": target.target,
        "image_backend": target.image_backend,
        "image_strategy": target.image_strategy,
        "variant": variant_name,
        "oci_dist": target.oci_dist,
        "base_image": base_image or target.base_image,
        "platforms": list(distribution.platforms),
        "context_dir": str(context_dir),
        "dockerfile": str(env_dockerfile),
        "bake_file": str(bake_file),
        "default_group": "all-dev",
        "available_groups": ["default", "all-dev", "all", "full", "full-all"],
        "component_targets": [spec.component_name for spec in component_specs],
        "docker_bake_command": f"docker buildx bake -f {bake_file} all-dev",
    }
