from __future__ import annotations

from pathlib import Path
import unittest

from feelpp.pkg.graph import build_plan, load_manifest, parse_skip_text


class GraphTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        repo_root = Path(__file__).resolve().parents[3]
        cls.manifest = load_manifest(repo_root / "packaging" / "manifest" / "components.toml")

    def test_default_order(self) -> None:
        plan = build_plan(self.manifest, dist="noble")
        self.assertEqual(
            [component.name for component in plan.components],
            ["feelpp", "feelpp-toolboxes", "feelpp-mor"],
        )

    def test_skip_text(self) -> None:
        selection = parse_skip_text("skip feelpp skip publish")
        self.assertIn("feelpp", selection.components)
        self.assertTrue(selection.publish)

    def test_requested_subset(self) -> None:
        plan = build_plan(
            self.manifest,
            dist="noble",
            requested_components=["feelpp-toolboxes", "feelpp-mor"],
            skipped_components={"feelpp-mor"},
        )
        self.assertEqual([component.name for component in plan.components], ["feelpp-toolboxes"])


if __name__ == "__main__":
    unittest.main()
