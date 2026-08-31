"""Check generated HTML links and assets inside the published pkgdown site.

No network access is used: absolute URLs under --base-url map to the same
build directory as relative URLs. External sites are outside this check.
"""

import argparse
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urljoin, urlsplit


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__(convert_charrefs=True)
        self.anchors = set()
        self.links = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get("id"):
            self.anchors.add(attrs["id"])
        if tag == "a" and attrs.get("name"):
            self.anchors.add(attrs["name"])
        attribute = "href" if tag in {"a", "link"} else "src"
        if tag in {"a", "link", "img", "script", "iframe", "source"}:
            value = attrs.get(attribute)
            if value:
                self.links.append((self.getpos()[0], value))


def check_site_links(directory, base_url="https://kfuku52.github.io/kfl1ou/"):
    root = Path(directory).resolve()
    if not (root / "index.html").is_file():
        raise ValueError(f"No generated index.html in {root}; build the site first.")
    base = urlsplit(base_url)
    if base.scheme not in {"http", "https"} or not base.netloc:
        raise ValueError("--base-url must be an absolute HTTP(S) site URL.")
    base_path = base.path.rstrip("/") + "/"
    base_url = f"{base.scheme}://{base.netloc}{base_path}"
    pages = {
        file.resolve(): Page(file.read_text(encoding="utf-8"))
        for file in root.rglob("*.html")
    }
    problems = []
    for source, page in sorted(pages.items()):
        source_name = source.relative_to(root).as_posix()
        source_url = urljoin(base_url, source_name)
        for line, href in page.links:
            target_url = urlsplit(urljoin(source_url, href))
            if (target_url.scheme not in {"http", "https"}
                    or target_url.netloc != base.netloc
                    or not target_url.path.startswith(base_path)):
                continue
            target = (root / unquote(target_url.path[len(base_path):])).resolve()
            if not target.is_relative_to(root):
                reason = "target escapes the site directory"
            else:
                if target.is_dir():
                    target = target / "index.html"
                if not target.is_file():
                    reason = "missing target"
                else:
                    fragment = unquote(target_url.fragment.split(":~:text=", 1)[0])
                    if target in pages and fragment and fragment not in pages[target].anchors:
                        reason = "missing anchor"
                    else:
                        continue
            problems.append(f"{source_name}:{line}: {href} ({reason})")
    return problems, len(pages)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", default="docs")
    parser.add_argument("--base-url", default="https://kfuku52.github.io/kfl1ou/")
    args = parser.parse_args()
    try:
        problems, pages = check_site_links(args.directory, args.base_url)
    except ValueError as error:
        parser.exit(1, f"{error}\n")
    if problems:
        parser.exit(1, "Broken internal links:\n" + "\n".join(problems) + "\n")
    print(f"Internal links and assets passed across {pages} HTML pages.")


if __name__ == "__main__":
    main()
