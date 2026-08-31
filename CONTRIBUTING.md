# Contributing to kfl1ou

Bug reports, reproducible examples, documentation fixes, and focused pull
requests are welcome.

## Development workflow

Open an issue first for changes that alter statistical behavior or the public
API. Work from `main` and keep each contribution focused. Add regression tests
for fixes and exercise new statistical behavior with deterministic examples.

During editing, use a focused test or a quick package check:

```sh
make bootstrap  # first checkout, or after dependency changes
make test FILTER=review-contracts
make check      # optional quick check; includes the complete test suite
```

Before submission, synchronize metadata and generated documentation, then run
the full package check. There is no need to run `make test` or `make check`
immediately before it; the full check already runs the tests.

```sh
git fetch origin
make bump        # synchronize version, UTC date, citation, NEWS, security line
make document    # after documentation/API edits
make metadata
make check-full  # includes vignettes and the PDF manual
```

For README or guide changes, also run `make doc-examples` to execute the actual
R examples with small simulation counts. Use `make docs` when changing the
site: it runs those examples, builds the site, and checks internal links and
assets. These checks share the existing documentation CI job.

Review and commit the intended changes, **including regenerated files**, then
run:

```sh
make generated
```

This final check regenerates documentation and requires all generated files to
match `HEAD`, including staged and previously untracked files. It is expected
to fail before intended generated changes are committed. See the
[development commands](https://github.com/kfuku52/kfl1ou/blob/main/tools/README.md)
for log retention, benchmarks and prerequisites.

Update `NEWS.md` when behavior visible to users changes. Every push must also
bump the version; `make bump` is idempotent relative to `origin/main` and
synchronizes release metadata. This includes dependency-only and Dependabot
updates. Review the generated NEWS entry. Use `make bump VERSION=3.1.0` when a
minor or major release is appropriate.

## Continuous integration

The matrix checks the supported R floor, current R on Linux, macOS and Windows,
coverage, generated files, dependency changes, and native code. Metadata
validation runs before expensive R checks. Changes to package code, tests,
build tools or dependencies run the complete matrix; metadata-only and
unrelated CI-action changes run the Linux package check. Manual and monthly
runs retain the full matrix, and releases dispatch omitted checks at the
verified tag before publishing. Scheduled extended checks cover R-devel,
oldrel, sanitizers, and Valgrind.

## Statistical changes

Describe the assumptions, criterion, parameterization, and numerical tolerance
of a statistical change. Whenever possible, include a small deterministic
simulation that would fail under the previous implementation. Avoid changing
defaults without documenting the compatibility and inferential consequences.

By contributing, you agree that your contribution is licensed under the
project's GPL-3.0-or-later license.
