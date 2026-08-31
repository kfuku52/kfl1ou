"""Regression fixtures for the documentation link checker (standard library only)."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from check_site_links import check_site_links


class SiteLinksTest(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory(prefix="kfl1ou-doc-links-")
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.write("index.html", '<h1 id="home">Home</h1>')

    def write(self, name, text):
        file = self.root / name
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_text(text, encoding="utf-8")

    def problems(self):
        return check_site_links(self.root)[0]

    def test_checks_relative_absolute_and_directory_links(self):
        self.write("articles/index.html", '<h1 id="guide">Guide</h1>')
        self.write("index.html", '''
            <a href="articles/#guide">relative</a>
            <a href="/kfl1ou/articles/index.html#guide">root relative</a>
            <a href="https://kfuku52.github.io/kfl1ou/articles/#guide">absolute</a>
            <a href="//kfuku52.github.io/kfl1ou/articles/#guide">protocol relative</a>
        ''')
        self.assertEqual(self.problems(), [])

    def test_catches_the_three_audited_missing_targets(self):
        self.write("index.html", '''
            <a href="LICENSE">license</a>
            <a href="NOTICE">notice</a>
            <a href="tools/README.html">development</a>
        ''')
        problems = self.problems()
        self.assertEqual(len(problems), 3)
        self.assertTrue(all("missing target" in p for p in problems))

    def test_catches_missing_fragments_and_accepts_encoded_anchors(self):
        self.write("index.html", '''
            <h1 id="with space">Heading</h1><a name="legacy"></a>
            <a href="#with%20space">encoded</a><a href="#legacy">legacy</a>
            <a href="#:~:text=Heading">text fragment</a><a href="#missing">bad</a>
        ''')
        self.assertEqual(len(self.problems()), 1)
        self.assertIn("#missing (missing anchor)", self.problems()[0])

    def test_checks_assets_and_ignores_queries(self):
        self.write("style.css", "body {}")
        self.write("image one.svg", "<svg/>")
        self.write("index.html", '''
            <link href="style.css?v=1&amp;theme=dark" rel="stylesheet">
            <img src="image%20one.svg"><script src="missing.js"></script>
        ''')
        self.assertEqual(len(self.problems()), 1)
        self.assertIn("missing.js", self.problems()[0])

    def test_external_targets_do_not_require_network_access(self):
        self.write("index.html", '''
            <a href="https://example.invalid/missing">external</a>
            <a href="mailto:example@example.invalid">email</a>
            <img src="data:image/png;base64,AA==">
            <a href="https://kfuku52.github.io/another-project/">another site</a>
        ''')
        self.assertEqual(self.problems(), [])

    def test_rejects_encoded_paths_that_escape_the_site(self):
        self.write("index.html", '<a href="%2e%2e/outside.html">escape</a>')
        self.assertIn("escapes the site directory", self.problems()[0])

    def test_empty_build_is_not_a_success(self):
        (self.root / "index.html").unlink()
        with self.assertRaisesRegex(ValueError, "build the site first"):
            self.problems()


if __name__ == "__main__":
    unittest.main()
