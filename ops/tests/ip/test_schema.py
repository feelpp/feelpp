from __future__ import annotations

import pytest

from feelpp.ops.ip.schema import validate_public_metadata

from .support import minimal_public_metadata


def test_valid_public_metadata_schema() -> None:
    validate_public_metadata(minimal_public_metadata())


def test_schema_rejects_private_metadata_keys() -> None:
    payload = minimal_public_metadata()
    payload["ownership_shares"] = []

    with pytest.raises(ValueError, match="Private metadata key"):
        validate_public_metadata(payload)


def test_schema_requires_repository_identity() -> None:
    payload = minimal_public_metadata()
    del payload["repository"]["default_branch"]

    with pytest.raises(ValueError, match="default_branch"):
        validate_public_metadata(payload)
