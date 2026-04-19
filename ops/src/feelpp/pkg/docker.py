from __future__ import annotations

from pathlib import Path
import hashlib
import json
import os
import re
import shlex
import shutil
from urllib.request import urlopen

from .apt_keys import (
    default_apt_keyring_keys,
    default_apt_signing_key,
    export_apt_public_keyring,
)
from .config import PackagingContext
from .shell import run, run_capture
from .workspace import ensure_workspace


DEFAULT_RUNTIME_PACKAGES = [
    "libfeelpp1",
    "feelpp-data",
    "feelpp-tools",
    "feelpp-quickstart",
    "python3-feelpp",
]
DEFAULT_APT_KEY_FILENAME = "feelpp-archive-keyring.gpg"
DEFAULT_REPO_CACHE_TOKEN_ENV = "FEELPP_PKG_IMAGE_REPO_TOKEN"
DEFAULT_OCI_REGISTRY = "ghcr.io"
DEFAULT_OCI_REPOSITORY = "feelpp/feelpp"
OCI_IMAGE_VERSION_LABEL = "org.opencontainers.image.version"
FEELPP_PACKAGE_VERSION_LABEL = "io.feelpp.package-version"
FEELPP_IMAGE_TAG_LABEL = "io.feelpp.image-tag"


def normalize_package_version_tag(package_version: str) -> str:
    normalized = package_version.split(":", 1)[-1]
    if "-" in normalized:
        normalized = normalized.rsplit("-", 1)[0]
    normalized = normalized.replace("~", "-")
    normalized = re.sub(r"[^a-zA-Z0-9_.-]+", "-", normalized).strip("-")
    if normalized and normalized[0].isdigit():
        normalized = f"v{normalized}"
    return normalized or "latest"


def runtime_image_tag_component(context: PackagingContext, *, feelpp_version: str) -> str:
    return f"{context.dist}-{normalize_package_version_tag(feelpp_version)}"


def default_runtime_image_tag(context: PackagingContext, *, feelpp_version: str) -> str:
    return f"feelpp:{runtime_image_tag_component(context, feelpp_version=feelpp_version)}"


def safe_workdir_name(raw: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_.-]+", "-", raw).strip("-")


def apt_release_metadata_url(context: PackagingContext) -> str:
    return f"http://apt.feelpp.org/{context.flavor}/{context.dist}/dists/{context.dist}/Release"


def resolve_repo_cache_token(
    context: PackagingContext,
    *,
    fallback: str,
    repo_cache_token: str | None = None,
) -> str:
    explicit_token = repo_cache_token or os.getenv(DEFAULT_REPO_CACHE_TOKEN_ENV)
    if explicit_token:
        return explicit_token

    try:
        with urlopen(apt_release_metadata_url(context), timeout=10) as response:
            release_bytes = response.read()
    except OSError:
        return fallback

    return hashlib.sha256(release_bytes).hexdigest()[:16]


def stage_runtime_apt_key(
    image_dir: Path,
    *,
    apt_key_file: str | None = None,
    gpg_key: str | None = None,
) -> Path:
    staged_key = image_dir / DEFAULT_APT_KEY_FILENAME
    if apt_key_file:
        source = Path(apt_key_file).expanduser().resolve()
        shutil.copyfile(source, staged_key)
        staged_key.chmod(0o644)
        return staged_key

    export_apt_public_keyring(
        staged_key,
        key_ids=default_apt_keyring_keys(primary_key=gpg_key),
    )
    staged_key.chmod(0o644)
    return staged_key


def render_runtime_dockerfile(
    context: PackagingContext,
    *,
    feelpp_version: str,
    packages: list[str],
    base_image: str,
    repo_cache_token_arg: str = "FEELPP_APT_REPO_TOKEN",
    apt_key_filename: str = DEFAULT_APT_KEY_FILENAME,
) -> str:
    quoted_packages = " ".join(shlex.quote(package) for package in packages)
    preference_packages = " ".join(packages)
    repo_url = f"http://apt.feelpp.org/{context.flavor}/{context.dist}"
    normalized_version = normalize_package_version_tag(feelpp_version)
    image_tag = runtime_image_tag_component(context, feelpp_version=feelpp_version)

    return f"""FROM {base_image}

ENV DEBIAN_FRONTEND=noninteractive
LABEL {OCI_IMAGE_VERSION_LABEL}="{normalized_version}" \\
      {FEELPP_PACKAGE_VERSION_LABEL}="{feelpp_version}" \\
      {FEELPP_IMAGE_TAG_LABEL}="{image_tag}"

RUN apt-get update \\
    && apt-get install -y --no-install-recommends ca-certificates \\
    && rm -rf /var/lib/apt/lists/*

RUN install -d -m 0755 /etc/apt/keyrings
COPY {apt_key_filename} /etc/apt/keyrings/feelpp.gpg
RUN chmod 0644 /etc/apt/keyrings/feelpp.gpg \\
    && printf '%s\\n' 'deb [signed-by=/etc/apt/keyrings/feelpp.gpg] {repo_url} {context.dist} {context.channel}' > /etc/apt/sources.list.d/feelpp.list

RUN cat >/etc/apt/preferences.d/feelpp <<'EOF'
Package: {preference_packages}
Pin: version {feelpp_version}
Pin-Priority: 1001
EOF

ARG {repo_cache_token_arg}=unknown
RUN apt-get update \\
    && printf '%s\\n' "${{{repo_cache_token_arg}}}" >/usr/local/share/feelpp-apt-repo-token \\
    && apt-get install -y --no-install-recommends {quoted_packages} \\
    && rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]
"""


def image_ref_tag(image_ref: str) -> str:
    tail = image_ref.rsplit("/", 1)[-1]
    if ":" in tail:
        return tail.rsplit(":", 1)[1]
    return "latest"


def inspect_docker_label(image_ref: str, label: str) -> str | None:
    value = run_capture(
        [
            "docker",
            "image",
            "inspect",
            "--format",
            f"{{{{index .Config.Labels {label!r}}}}}",
            image_ref,
        ],
        check=False,
    )
    value = value.strip()
    if not value or value == "<no value>":
        return None
    return value


def resolve_publish_tag(source_ref: str, *, tag: str | None = None) -> str:
    if tag:
        return tag
    labeled_tag = inspect_docker_label(source_ref, FEELPP_IMAGE_TAG_LABEL)
    if labeled_tag:
        return labeled_tag
    return image_ref_tag(source_ref)


def default_oci_registry() -> str:
    return os.getenv("FEELPP_PKG_OCI_REGISTRY") or DEFAULT_OCI_REGISTRY


def default_oci_repository() -> str:
    return os.getenv("FEELPP_PKG_OCI_REPOSITORY") or DEFAULT_OCI_REPOSITORY


def default_publish_ref(
    source_ref: str,
    *,
    target_ref: str | None = None,
    registry: str | None = None,
    repository: str | None = None,
    tag: str | None = None,
) -> str:
    if target_ref:
        return target_ref
    resolved_registry = registry or default_oci_registry()
    resolved_repository = repository or default_oci_repository()
    resolved_tag = resolve_publish_tag(source_ref, tag=tag)
    return f"{resolved_registry}/{resolved_repository}:{resolved_tag}"


def build_runtime_image(
    context: PackagingContext,
    *,
    feelpp_version: str,
    tag: str | None = None,
    packages: list[str] | None = None,
    base_image: str | None = None,
    apt_key_file: str | None = None,
    gpg_key: str | None = None,
    repo_cache_token: str | None = None,
    no_cache: bool = False,
    dry_run: bool = False,
) -> dict[str, object]:
    ensure_workspace(context)
    resolved_packages = packages or list(DEFAULT_RUNTIME_PACKAGES)
    resolved_base_image = base_image or f"{context.flavor}:{context.dist}"
    resolved_tag = tag or default_runtime_image_tag(context, feelpp_version=feelpp_version)
    resolved_repo_cache_token = resolve_repo_cache_token(
        context,
        fallback=feelpp_version,
        repo_cache_token=repo_cache_token,
    )

    image_dir = context.job_root / "images" / safe_workdir_name(resolved_tag)
    image_dir.mkdir(parents=True, exist_ok=True)
    staged_key_path = stage_runtime_apt_key(
        image_dir,
        apt_key_file=apt_key_file,
        gpg_key=gpg_key,
    )
    dockerfile_path = image_dir / "Dockerfile"
    dockerfile_path.write_text(
        render_runtime_dockerfile(
            context,
            feelpp_version=feelpp_version,
            packages=resolved_packages,
            base_image=resolved_base_image,
            apt_key_filename=staged_key_path.name,
        ),
        encoding="utf-8",
    )

    command = ["docker", "build"]
    if no_cache:
        command.append("--no-cache")
    command.extend(
        [
            "--build-arg",
            f"FEELPP_APT_REPO_TOKEN={resolved_repo_cache_token}",
            "-t",
            resolved_tag,
            "-f",
            str(dockerfile_path),
            str(image_dir),
        ]
    )
    run(command, cwd=context.repo_root, env=context.shell_env(), dry_run=dry_run)

    result = {
        "tag": resolved_tag,
        "version_tag": runtime_image_tag_component(context, feelpp_version=feelpp_version),
        "base_image": resolved_base_image,
        "feelpp_version": feelpp_version,
        "packages": resolved_packages,
        "repo_cache_token": resolved_repo_cache_token,
        "repo_release_url": apt_release_metadata_url(context),
        "apt_key": str(staged_key_path),
        "dockerfile": str(dockerfile_path),
        "build_context": str(image_dir),
    }
    print(json.dumps(result, indent=2))
    return result


def publish_docker_image(
    context: PackagingContext,
    *,
    source_ref: str,
    target_ref: str | None = None,
    registry: str | None = None,
    repository: str | None = None,
    tag: str | None = None,
    dry_run: bool = False,
) -> dict[str, str]:
    ensure_workspace(context)
    resolved_target_ref = default_publish_ref(
        source_ref,
        target_ref=target_ref,
        registry=registry,
        repository=repository,
        tag=tag,
    )
    run(
        ["docker", "tag", source_ref, resolved_target_ref],
        cwd=context.repo_root,
        env=context.shell_env(),
        dry_run=dry_run,
    )
    run(
        ["docker", "push", resolved_target_ref],
        cwd=context.repo_root,
        env=context.shell_env(),
        dry_run=dry_run,
    )
    result = {
        "source_ref": source_ref,
        "target_ref": resolved_target_ref,
    }
    print(json.dumps(result, indent=2))
    return result
