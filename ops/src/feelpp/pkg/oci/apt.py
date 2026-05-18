from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import shutil

from jinja2 import Environment, FileSystemLoader, StrictUndefined

from .catalog import ImageTarget, list_image_targets
from .bake_permissions import fs_read_allow_flags, required_fs_read_paths_from_bake_payload
from .cmake_presets import resolve_cmake_preset
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
from .platforms import resolve_target_platforms
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

COMPONENT_REQUEST_ALIASES = {
    "env": "env",
    "feelpp": "feelpp",
    "toolboxes": "toolboxes",
    "feelpp-toolboxes": "toolboxes",
    "mor": "mor",
    "feelpp-mor": "mor",
    "full": "full",
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


def _normalize_requested_component(component: str | None) -> str | None:
    if component in {None, ""}:
        return None
    normalized = COMPONENT_REQUEST_ALIASES.get(str(component).strip().lower())
    if not normalized:
        known = ", ".join(sorted(COMPONENT_REQUEST_ALIASES))
        raise ValueError(f"Unsupported image component: {component}. Known values: {known}")
    return normalized


def _component_context_dir_name(component_name: str) -> str:
    return {
        "feelpp": "feelpp",
        "feelpp-toolboxes": "feelpp-toolboxes",
        "feelpp-mor": "feelpp-mor",
    }[component_name]


def _component_image_names(component_name: str) -> tuple[str, str]:
    if component_name == "feelpp":
        return "feelpp", "feelpp"
    if component_name == "feelpp-toolboxes":
        return "feelpp-toolboxes", "feelpp-toolboxes"
    return "feelpp-mor", "feelpp-mor"


def _selected_groups(component: str | None) -> list[str]:
    if component in {None, "env"}:
        return ["default"]
    if component == "feelpp":
        return ["feelpp-all"]
    if component == "toolboxes":
        return ["toolboxes-all"]
    if component == "mor":
        return ["mor-all"]
    if component == "full":
        return ["full-all"]
    raise ValueError(f"Unsupported image component selection: {component}")


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
    component: str | None = None,
    from_image: str | None = None,
    platform_overrides: list[str] | None = None,
) -> dict[str, object]:
    metadata = load_docker_metadata(workspace.repo_root)
    distribution = resolve_distribution_version(metadata, family=target.flavor, name=target.dist)
    component_specs = _component_targets_for_row(workspace.repo_root, target)
    requested_component = _normalize_requested_component(component)
    resolved_platforms = resolve_target_platforms(distribution.platforms, platform_overrides)

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
            "platforms": resolved_platforms,
        }
    }
    image_refs: dict[str, str] = {
        "feelpp-env": _environment_image_ref(target, registry=registry, namespace=namespace),
    }

    previous_runtime_target = "feelpp-env"
    previous_context_alias = "feelpp_env_image"
    for component_spec in component_specs:
        context_dir_name = _component_context_dir_name(component_spec.component_name)
        builder_tag = _builder_tag(target, branch=workspace.branch)
        runtime_tag = _runtime_tag(target, branch=workspace.branch)
        builder_image_name, runtime_image_name = _component_image_names(component_spec.component_name)

        builder_ref = _component_image_ref(
            builder_image_name,
            builder_tag,
            registry=registry,
            namespace=namespace,
        )
        runtime_ref = _component_image_ref(
            runtime_image_name,
            runtime_tag,
            registry=registry,
            namespace=namespace,
        )
        image_refs[f"{component_spec.bake_target}-dev"] = builder_ref
        image_refs[component_spec.bake_target] = runtime_ref

        contexts = {"feelpp_source": str(workspace.repo_root)}
        from_image_arg = previous_context_alias
        if requested_component == component_spec.bake_target and from_image:
            from_image_arg = from_image
        else:
            contexts[previous_context_alias] = f"target:{previous_runtime_target}"
        runtime_contexts = dict(contexts)
        runtime_contexts["build_output"] = f"target:{component_spec.bake_target}"
        builder_args = {
            "FROM_IMAGE": from_image_arg,
            "DESCRIPTION": f"{component_spec.description} (dev)",
            "BRANCH": workspace.branch,
            "CMAKE_PRESET": resolve_cmake_preset(target, component=component_spec.bake_target),
            "CXX": cxx,
            "CC": cc,
            "CMAKE_FLAGS": cmake_flags,
        }
        runtime_args = {
            "FROM_IMAGE": from_image_arg,
            "DESCRIPTION": component_spec.description,
            "BRANCH": workspace.branch,
            "CMAKE_PRESET": resolve_cmake_preset(target, component=component_spec.bake_target),
            "CXX": cxx,
            "CC": cc,
            "CMAKE_FLAGS": cmake_flags,
        }
        if component_spec.bake_target == "toolboxes":
            builder_args["RUN_CTEST"] = "1"
            runtime_args["RUN_CTEST"] = "1"

        targets[component_spec.bake_target] = {
            "context": context_dir_name,
            "dockerfile": "Dockerfile.multistage",
            "target": "builder",
            "contexts": contexts,
            "args": builder_args,
            "tags": [builder_ref],
            "platforms": resolved_platforms,
        }
        targets[component_spec.runtime_target] = {
            "context": context_dir_name,
            "dockerfile": "Dockerfile.multistage",
            "target": "runtime",
            "contexts": runtime_contexts,
            "args": runtime_args,
            "tags": [runtime_ref],
            "platforms": resolved_platforms,
        }
        previous_runtime_target = component_spec.runtime_target
        previous_context_alias = f"{component_spec.bake_target}_runtime_image"

    full_builder_ref = _component_image_ref(
        "feelpp",
        _full_builder_tag(target, branch=workspace.branch),
        registry=registry,
        namespace=namespace,
    )
    full_runtime_ref = _component_image_ref(
        "feelpp",
        _full_runtime_tag(target, branch=workspace.branch),
        registry=registry,
        namespace=namespace,
    )
    image_refs["feelpp-full-dev"] = full_builder_ref
    image_refs["feelpp-full"] = full_runtime_ref

    full_contexts = {"feelpp_source": str(workspace.repo_root)}
    full_from_image_arg = "feelpp_env_image"
    if requested_component == "full" and from_image:
        full_from_image_arg = from_image
    else:
        full_contexts["feelpp_env_image"] = "target:feelpp-env"
    full_runtime_contexts = dict(full_contexts)
    full_runtime_contexts["build_output"] = "target:feelpp-full"
    targets["feelpp-full"] = {
        "context": "feelpp",
        "dockerfile": "Dockerfile.full",
        "target": "builder",
        "contexts": full_contexts,
        "args": {
            "FROM_IMAGE": full_from_image_arg,
            "DESCRIPTION": "Feel++ Full Stack (dev)",
            "BRANCH": workspace.branch,
            "CMAKE_PRESET": resolve_cmake_preset(target, component="full"),
            "CXX": cxx,
            "CC": cc,
            "CMAKE_FLAGS": cmake_flags,
        },
        "tags": [full_builder_ref],
        "platforms": resolved_platforms,
    }
    targets["feelpp-full-runtime"] = {
        "context": "feelpp",
        "dockerfile": "Dockerfile.full",
        "target": "runtime",
        "contexts": full_runtime_contexts,
        "args": {
            "FROM_IMAGE": full_from_image_arg,
            "DESCRIPTION": "Feel++ Full Stack",
            "BRANCH": workspace.branch,
            "CMAKE_PRESET": resolve_cmake_preset(target, component="full"),
            "CXX": cxx,
            "CC": cc,
            "CMAKE_FLAGS": cmake_flags,
        },
        "tags": [full_runtime_ref],
        "platforms": resolved_platforms,
    }

    bake_payload = {
        "group": {
            "default": {"targets": ["feelpp-env"]},
            "env": {"targets": ["feelpp-env"]},
            "feelpp-dev": {"targets": ["feelpp"]},
            "feelpp-runtime": {"targets": ["feelpp-runtime"]},
            "feelpp-all": {"targets": ["feelpp", "feelpp-runtime"]},
            "toolboxes-dev": {"targets": ["toolboxes"]},
            "toolboxes-runtime": {"targets": ["toolboxes-runtime"]},
            "toolboxes-all": {"targets": ["toolboxes", "toolboxes-runtime"]},
            "mor-dev": {"targets": ["mor"]},
            "mor-runtime": {"targets": ["mor-runtime"]},
            "mor-all": {"targets": ["mor", "mor-runtime"]},
            "all-dev": {"targets": [spec.bake_target for spec in component_specs]},
            "all": {"targets": [spec.runtime_target for spec in component_specs]},
            "full-dev": {"targets": ["feelpp-full"]},
            "full-runtime": {"targets": ["feelpp-full-runtime"]},
            "full": {"targets": ["feelpp-full"]},
            "full-all": {"targets": ["feelpp-full", "feelpp-full-runtime"]},
        },
        "target": targets,
    }
    bake_file = context_dir / "docker-bake.json"
    bake_file.write_text(json.dumps(bake_payload, indent=2) + "\n", encoding="utf-8")
    recommended_groups = _selected_groups(requested_component)
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
        "variant": variant_name,
        "oci_dist": target.oci_dist,
        "base_image": base_image or target.base_image,
        "selected_component": requested_component or "env",
        "from_image": from_image or "",
        "platforms": resolved_platforms,
        "context_dir": str(context_dir),
        "dockerfile": str(env_dockerfile),
        "bake_file": str(bake_file),
        "default_group": "default",
        "available_groups": [
            "default",
            "env",
            "feelpp-dev",
            "feelpp-runtime",
            "feelpp-all",
            "toolboxes-dev",
            "toolboxes-runtime",
            "toolboxes-all",
            "mor-dev",
            "mor-runtime",
            "mor-all",
            "all-dev",
            "all",
            "full-dev",
            "full-runtime",
            "full",
            "full-all",
        ],
        "recommended_groups": recommended_groups,
        "fs_read_paths": fs_read_paths,
        "component_targets": [spec.component_name for spec in component_specs],
        "image_refs": image_refs,
        "docker_bake_command": " ".join(docker_bake_command),
    }
