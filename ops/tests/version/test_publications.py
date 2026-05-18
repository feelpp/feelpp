from __future__ import annotations

from unittest import mock
import io
import json
import tempfile
import unittest

from feelpp.ops.version.cli import build_parser as build_version_parser, main as version_main
from feelpp.ops.version.publications import HalPublication, HalPublicationService

from ops.tests.version.support import VersionRepoMixin


class PublicationsTests(VersionRepoMixin, unittest.TestCase):
    def test_hal_publication_service_uses_feel_and_cemosis_collections(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            service = HalPublicationService(repo_root)
            payload = {
                "response": {
                    "docs": [
                        {
                            "docid": "1",
                            "title_s": ["A Feel++ paper"],
                            "uri_s": "https://hal.science/hal-1",
                            "producedDate_tdate": "2026-01-28T00:00:00Z",
                            "producedDateY_i": 2026,
                            "docType_s": "UNDEFINED",
                            "authFullName_s": ["Alice", "Bob", "Carol", "Dave"],
                            "journalTitle_s": "HAL Preprint",
                        }
                    ]
                }
            }

            class FakeResponse:
                def __enter__(self):
                    return self

                def __exit__(self, exc_type, exc, tb):
                    return False

                def read(self):
                    return json.dumps(payload).encode("utf-8")

            with mock.patch("feelpp.ops.version.publications.urlopen", return_value=FakeResponse()) as urlopen_mock:
                publications = service.fetch()

            request = urlopen_mock.call_args.args[0]
            self.assertIn("collCode_s%3A%28FEEL+OR+CEMOSIS%29", request.full_url)
            self.assertEqual(len(publications), 1)
            self.assertEqual(publications[0].title, "A Feel++ paper")
            markdown = service.format_markdown(publications)
            self.assertIn("## Recent Publications using Feel++", markdown)
            self.assertIn("Alice, Bob, Carol, et al.", markdown)
            self.assertIn("Preprint", markdown)

    def test_cli_parses_release_publication_flags(self) -> None:
        parser = build_version_parser()
        args = parser.parse_args(
            [
                "release",
                "0.111.0-preview.13",
                "--publications-rows",
                "7",
                "--publications-since",
                "2026-01-01",
            ]
        )
        self.assertEqual(args.publications_rows, 7)
        self.assertEqual(args.publications_since, "2026-01-01")

    def test_cli_parses_publications_pretty_flag(self) -> None:
        parser = build_version_parser()
        args = parser.parse_args(["publications", "--pretty", "--collection", "FEEL"])
        self.assertTrue(args.pretty)
        self.assertEqual(args.collection, ["FEEL"])

    def test_publications_command_pretty_prints_markdown(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            stdout = io.StringIO()
            publication = HalPublication(
                docid="1",
                title="A Feel++ paper",
                url="https://hal.science/hal-1",
                produced_at="2026-01-28T00:00:00Z",
                year=2026,
                doc_type="UNDEFINED",
                venue="HAL Preprint",
                authors=("Alice", "Bob"),
                doi=None,
            )

            with mock.patch("feelpp.ops.version.cli.HalPublicationService") as service_cls:
                service_cls.return_value.fetch.return_value = [publication]
                service_cls.return_value.format_markdown.return_value = "## Recent Publications using Feel++\n\n- [A Feel++ paper](https://hal.science/hal-1)"
                with mock.patch("sys.stdout", stdout):
                    rc = version_main(["--repo-root", str(repo_root), "publications", "--pretty"])

            self.assertEqual(rc, 0)
            self.assertIn("## Recent Publications using Feel++", stdout.getvalue())

    def test_release_command_passes_publication_overrides(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            repo_root = self.make_repo(tmpdir)
            stdout = io.StringIO()

            class FakePlan:
                title = "v0.111.0-preview.13"
                tag = "v0.111.0-preview.13"
                release_kind = "upstream"
                branch = "develop"
                channel = "latest"
                head_sha = "abc123"
                previous_tag = "v0.111.0-preview.12"
                prerelease = True
                dry_run = True
                package_notes = "## Packages"
                generated_notes_preview = "## What's Changed"

            with mock.patch("feelpp.ops.version.cli.ReleaseService") as service_cls:
                service_cls.return_value.execute_release.return_value = FakePlan()
                with mock.patch("sys.stdout", stdout):
                    rc = version_main(
                        [
                            "--repo-root",
                            str(repo_root),
                            "release",
                            "0.111.0-preview.13",
                            "--publications-rows",
                            "7",
                            "--publications-since",
                            "2026-01-01",
                            "--pretty",
                        ]
                    )

            self.assertEqual(rc, 0)
            service_cls.return_value.execute_release.assert_called_once_with(
                "0.111.0-preview.13",
                dry_run=False,
                dists=None,
                publication_rows=7,
                publication_since="2026-01-01",
            )


if __name__ == "__main__":
    unittest.main()

