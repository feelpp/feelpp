from __future__ import annotations

import importlib

from feelpp.ops.common import (
    PREFERRED_DEV_CLI_NAME,
    LEGACY_PACKAGE_CLI_NAMES,
    LEGACY_VERSION_CLI_NAMES,
    PACKAGE_CLI_ALIAS_COMMANDS,
    PREFERRED_PACKAGE_CLI_NAME,
    PREFERRED_VERSION_CLI_NAME,
)
from feelpp.ops.dev.cli import build_parser as build_dev_parser
from feelpp.pkg.cli import build_parser
from feelpp.ops.version.cli import build_parser as build_version_parser


def test_pkg_cli_uses_shared_preferred_name() -> None:
    assert build_parser().prog == PREFERRED_PACKAGE_CLI_NAME
    assert LEGACY_PACKAGE_CLI_NAMES == ("feelpp-pkg",)
    assert PACKAGE_CLI_ALIAS_COMMANDS == {"fpp-spack": ("spack",)}


def test_pkg_cli_supports_explicit_backend_groups() -> None:
    parser = build_parser()
    debian_args = parser.parse_args(["debian", "job", "init", "--dist", "noble"])
    spack_args = parser.parse_args(["spack", "env", "list"])
    assert debian_args.command == "debian"
    assert spack_args.command == "spack"


def test_version_cli_uses_shared_preferred_name() -> None:
    assert build_version_parser().prog == PREFERRED_VERSION_CLI_NAME
    assert LEGACY_VERSION_CLI_NAMES == ("feelpp-version",)


def test_dev_cli_uses_shared_preferred_name() -> None:
    assert build_dev_parser().prog == PREFERRED_DEV_CLI_NAME


def test_ops_common_namespace_is_importable() -> None:
    module = importlib.import_module("feelpp.ops.common")
    assert module.PREFERRED_PACKAGE_CLI_NAME == "fpp-pkg"
    assert module.PACKAGE_CLI_ALIAS_COMMANDS == {"fpp-spack": ("spack",)}
    assert module.PREFERRED_DEV_CLI_NAME == "fpp-dev"


def test_pkg_backend_subpackages_are_importable() -> None:
    debian_backend = importlib.import_module("feelpp.pkg.backends.debian")
    spack_backend = importlib.import_module("feelpp.pkg.backends.spack")
    core_module = importlib.import_module("feelpp.pkg.core")
    assert debian_backend.__name__ == "feelpp.pkg.backends.debian"
    assert spack_backend.__name__ == "feelpp.pkg.backends.spack"
    assert core_module.__name__ == "feelpp.pkg.core"


def test_simulate_namespace_is_importable() -> None:
    module = importlib.import_module("feelpp.simulate")
    assert module.__name__ == "feelpp.simulate"
