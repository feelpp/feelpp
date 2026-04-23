from __future__ import annotations

import importlib

from feelpp.ops.common import (
    LEGACY_PACKAGE_CLI_NAMES,
    LEGACY_VERSION_CLI_NAMES,
    PREFERRED_IP_CLI_NAME,
    PREFERRED_PACKAGE_CLI_NAME,
    PREFERRED_VERSION_CLI_NAME,
)
from feelpp.ops.ip.cli import build_parser as build_ip_parser
from feelpp.pkg.cli import build_parser
from feelpp.ops.version.cli import build_parser as build_version_parser


def test_pkg_cli_uses_shared_preferred_name() -> None:
    assert build_parser().prog == PREFERRED_PACKAGE_CLI_NAME
    assert LEGACY_PACKAGE_CLI_NAMES == ("feelpp-pkg",)


def test_version_cli_uses_shared_preferred_name() -> None:
    assert build_version_parser().prog == PREFERRED_VERSION_CLI_NAME
    assert LEGACY_VERSION_CLI_NAMES == ("feelpp-version",)


def test_ip_cli_uses_shared_preferred_name() -> None:
    assert build_ip_parser().prog == PREFERRED_IP_CLI_NAME


def test_ops_common_namespace_is_importable() -> None:
    module = importlib.import_module("feelpp.ops.common")
    assert module.PREFERRED_PACKAGE_CLI_NAME == "fpp-pkg"
    assert module.PREFERRED_IP_CLI_NAME == "fpp-ip"


def test_simulate_namespace_is_importable() -> None:
    module = importlib.import_module("feelpp.simulate")
    assert module.__name__ == "feelpp.simulate"
