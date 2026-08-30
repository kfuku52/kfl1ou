# Contributing to kfl1ou

Bug reports, reproducible examples, documentation fixes, and focused pull
requests are welcome.

## Before opening a pull request

1. Open an issue first for changes that alter statistical behavior or the
   public API.
2. Work from the `main` branch and keep each pull request focused on one
   coherent change.
3. Add regression tests for fixes and tests that exercise every new branch of
   statistical code.
4. Use a focused test while editing, then run one complete package check:

   ```sh
   make bootstrap  # first checkout, or after dependency changes
   make test FILTER=review-contracts  # only the affected test file(s)
   make check       # all tests, without rebuilding vignettes or the PDF manual
   make bump        # synchronize version, UTC date, citation, NEWS, security line
   make document    # update generated files after documentation/API edits
   make generated
   make metadata
   make check-full  # before submitting
   ```

5. Update `NEWS.md` when behavior visible to users changes. Every push must also
   bump the version; `make bump` updates the related metadata together and is
   idempotent relative to `origin/main`. This also applies to dependency-only
   updates, including Dependabot updates. Review the generated NEWS entry.
   Use `make bump VERSION=3.1.0` when a minor or major release is appropriate.

The continuous-integration matrix checks the supported R floor, current R on
Linux, macOS and Windows, coverage, generated files, dependency changes, and
native code. Metadata validation runs before expensive R checks. Changes to
package code, tests, build tools or dependencies run the complete matrix;
metadata-only and unrelated CI-action changes run the Linux package check.
Manual and monthly runs retain the full matrix, and releases dispatch any
omitted checks at the verified tag before publishing. Scheduled extended
checks cover R-devel, oldrel, sanitizers, and Valgrind.

`make test` without `FILTER` runs the complete test suite. It is useful for
debugging, but is not required immediately before `make check`, which already
runs those tests. See [development commands](tools/README.md) for benchmarks
and metadata maintenance.

## Statistical changes

Describe the assumptions, criterion, parameterization, and numerical tolerance
of a statistical change. Whenever possible, include a small deterministic
simulation that would fail under the previous implementation. Avoid changing
defaults without documenting the compatibility and inferential consequences.

By contributing, you agree that your contribution is licensed under the
project's GPL-3.0-or-later license.
