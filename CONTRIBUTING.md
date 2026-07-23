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
4. Run the package tests and check generated files:

   ```sh
   Rscript -e 'devtools::test()'
   Rscript -e 'Rcpp::compileAttributes("."); roxygen2::roxygenise(".")'
   git diff --exit-code
   R CMD build .
   R CMD check --as-cran kfl1ou_*.tar.gz
   ```

5. Update `NEWS.md` when behavior visible to users changes.

The continuous-integration matrix checks the supported R floor, current R on
Linux, macOS and Windows, coverage, generated files, dependency changes, and
native code. Scheduled extended checks cover R-devel, oldrel, sanitizers, and
Valgrind.

## Statistical changes

Describe the assumptions, criterion, parameterization, and numerical tolerance
of a statistical change. Whenever possible, include a small deterministic
simulation that would fail under the previous implementation. Avoid changing
defaults without documenting the compatibility and inferential consequences.

By contributing, you agree that your contribution is licensed under the
project's GPL-3.0-or-later license.
