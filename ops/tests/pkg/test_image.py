from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest import mock

from feelpp.pkg.config import PackagingContext
from feelpp.pkg.image import (
    DEFAULT_APT_KEY_FILENAME,
    DEFAULT_RUNTIME_PACKAGES,
    apt_release_metadata_url,
    build_runtime_image,
    default_apt_signing_key,
    default_runtime_image_tag,
    render_runtime_dockerfile,
    resolve_repo_cache_token,
    stage_runtime_apt_key,
)


class ImageTests(unittest.TestCase):
    def make_context(self, tmpdir: str) -> PackagingContext:
        repo_root = Path(tmpdir) / "repo"
        (repo_root / "packaging" / "pbuilder" / "hooks").mkdir(parents=True)
        (repo_root / "packaging" / "pbuilder" / "pbuilderrc").write_text("", encoding="utf-8")
        job_root = repo_root / "job"
        return PackagingContext.create(
            repo_root=repo_root,
            dist="noble",
            flavor="ubuntu",
            branch="develop",
            channel="latest",
            job_id="test-job",
            job_root=job_root,
        )

    def test_default_runtime_image_tag(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            self.assertEqual(
                default_runtime_image_tag(context, feelpp_version="0.111.0~preview.13-1"),
                "feelpp:noble-v0.111.0-preview.13",
            )

    def test_default_apt_signing_key_prefers_env(self) -> None:
        with mock.patch.dict("os.environ", {"GPG_KEY": "TESTKEY"}, clear=False):
            self.assertEqual(default_apt_signing_key(), "TESTKEY")

    def test_render_runtime_dockerfile_pins_requested_version(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            dockerfile = render_runtime_dockerfile(
                context,
                feelpp_version="0.111.0~preview.13-1",
                packages=["libfeelpp1", "python3-feelpp"],
                base_image="ubuntu:noble",
            )

            self.assertIn("FROM ubuntu:noble", dockerfile)
            self.assertIn(f"COPY {DEFAULT_APT_KEY_FILENAME} /etc/apt/keyrings/feelpp.gpg", dockerfile)
            self.assertIn("deb [signed-by=/etc/apt/keyrings/feelpp.gpg] http://apt.feelpp.org/ubuntu/noble noble latest", dockerfile)
            self.assertIn("Pin: version 0.111.0~preview.13-1", dockerfile)
            self.assertIn("ARG FEELPP_APT_REPO_TOKEN=unknown", dockerfile)
            self.assertIn("printf '%s\\n' \"${FEELPP_APT_REPO_TOKEN}\" >/usr/local/share/feelpp-apt-repo-token", dockerfile)
            self.assertIn("apt-get install -y --no-install-recommends libfeelpp1 python3-feelpp", dockerfile)

    def test_apt_release_metadata_url(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)
            self.assertEqual(
                apt_release_metadata_url(context),
                "http://apt.feelpp.org/ubuntu/noble/dists/noble/Release",
            )

    def test_resolve_repo_cache_token_hashes_release_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)

            class FakeResponse:
                def __enter__(self):
                    return self

                def __exit__(self, exc_type, exc, tb):
                    return False

                def read(self) -> bytes:
                    return b"Date: Wed, 1 Apr 2026 03:15:05 UTC\n"

            with mock.patch("feelpp.pkg.docker.urlopen", return_value=FakeResponse()):
                token = resolve_repo_cache_token(context, fallback="0.111.0~preview.13-1")

            self.assertEqual(token, "f1a08aa0d28724a4")

    def test_resolve_repo_cache_token_falls_back_when_release_unavailable(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)

            with mock.patch("feelpp.pkg.docker.urlopen", side_effect=OSError("offline")):
                token = resolve_repo_cache_token(context, fallback="0.111.0~preview.13-1")

            self.assertEqual(token, "0.111.0~preview.13-1")

    def test_stage_runtime_apt_key_copies_explicit_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            image_dir = Path(tmpdir) / "image"
            image_dir.mkdir()
            source_key = Path(tmpdir) / "feelpp.gpg"
            source_key.write_bytes(b"keydata")

            staged = stage_runtime_apt_key(image_dir, apt_key_file=str(source_key))

            self.assertEqual(staged.read_bytes(), b"keydata")
            self.assertEqual(staged.name, DEFAULT_APT_KEY_FILENAME)

    def test_build_runtime_image_writes_dockerfile_and_invokes_docker(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)

            source_key = Path(tmpdir) / "feelpp.gpg"
            source_key.write_bytes(b"keydata")
            with mock.patch("feelpp.pkg.docker.run") as run:
                with mock.patch("feelpp.pkg.docker.resolve_repo_cache_token", return_value="repo-release-token"):
                    result = build_runtime_image(
                        context,
                        feelpp_version="0.111.0~preview.13-1",
                        apt_key_file=str(source_key),
                        dry_run=True,
                    )

            self.assertEqual(result["packages"], DEFAULT_RUNTIME_PACKAGES)
            self.assertEqual(result["repo_cache_token"], "repo-release-token")
            dockerfile_path = Path(result["dockerfile"])
            self.assertTrue(dockerfile_path.is_file())
            self.assertIn("0.111.0~preview.13-1", dockerfile_path.read_text(encoding="utf-8"))
            self.assertTrue(Path(result["apt_key"]).is_file())

            run.assert_called_once()
            command = run.call_args.args[0]
            self.assertEqual(command[:2], ["docker", "build"])
            self.assertIn("--build-arg", command)
            self.assertIn("FEELPP_APT_REPO_TOKEN=repo-release-token", command)
            self.assertIn(str(dockerfile_path), command)

    def test_build_runtime_image_supports_no_cache(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            context = self.make_context(tmpdir)

            source_key = Path(tmpdir) / "feelpp.gpg"
            source_key.write_bytes(b"keydata")
            with mock.patch("feelpp.pkg.docker.run") as run:
                with mock.patch("feelpp.pkg.docker.resolve_repo_cache_token", return_value="repo-release-token"):
                    result = build_runtime_image(
                        context,
                        feelpp_version="0.111.0~preview.13-1",
                        apt_key_file=str(source_key),
                        no_cache=True,
                        dry_run=True,
                    )

            run.assert_called_once()
            command = run.call_args.args[0]
            self.assertEqual(command[:3], ["docker", "build", "--no-cache"])
            self.assertIn("FEELPP_APT_REPO_TOKEN=repo-release-token", command)


if __name__ == "__main__":
    unittest.main()
