# Development commands

Run these commands from the repository root. `make bootstrap` installs the
development dependencies; the package itself still has only its declared
runtime dependencies.
The full check also needs system TeX tools (including `makeindex`) and `qpdf`.
Use an HTML5-capable `tidy` when validating the generated HTML reference manual.
Maintenance fixtures use Python >= 3.9, Bash and `jq`, as available on hosted CI.
Documentation builds need Pandoc; `make bootstrap` also installs pkgdown.

| Command | Purpose |
|---|---|
| `make test FILTER=measurement-error` | Run test files matching a regular expression. |
| `make test` | Run the complete test suite without a package check. |
| `make check` | Build and check the package, including tests and examples, without rebuilding vignettes or the PDF manual. |
| `make check-full` | Build vignettes and run `R CMD check --as-cran`, including the PDF manual; requires TeX. |
| `make document` | Regenerate Rcpp exports and Rd documentation. |
| `make generated` | Regenerate and require generated files to match HEAD, including staged and untracked changes; run after committing reviewed output. |
| `make doc-examples` | Run the README quick start and guide examples, including normally skipped calibration blocks with two simulations. |
| `make docs` | Run documentation examples, build the pkgdown site and check internal links and assets. |
| `make bump` | Ensure a patch version above `origin/main` and synchronize all release metadata. |
| `make metadata` | Check metadata consistency and maintenance-script fixtures. |
| `make benchmark` | Install this checkout in a temporary library and run deterministic smoke benchmarks. |
| `KFL1OU_BENCHMARK_PROFILE=extended make benchmark` | Add the larger explicit-pruning workload. |

Use focused tests during implementation and `make check-full` before submission.
There is no need to run the full test suite immediately before a package check.
Set `KFL1OU_CHECK_OUTPUT` to a directory to retain build and check logs there;
otherwise temporary check files are removed at exit.
After changing documentation, run `make document`, review and commit its output,
and then use `make generated` to verify reproducibility. The check compares
against HEAD, not just the index, and detects new untracked generated pages.
Intended generated changes must therefore be committed before this final check.

`make doc-examples` installs this checkout in a temporary library, including
for socket workers, then executes the actual named R blocks in the guides with the
same search settings, reducing only the simulation count. This is an execution
check, not statistical calibration. `make docs` also builds the site and checks
relative links, absolute links under the site URL, anchors and local assets.
The site build reuses that installation. CI uses `--installed` only after
installing the checkout through `setup-r-dependencies`. External links are
outside the offline link check. To check an existing build, run
`python3 tools/check_site_links.py docs`; use `--base-url` for another deployment.

Keep the README quick-start markers and named guide chunks in sync with
`tools/check-doc-examples.R`. A missing chunk is an error rather than a silently
skipped example. `make metadata` runs the small link-checker and generated-file
fixtures; the pkgdown job runs the examples and checks its built artifact before
publishing, without adding another OS matrix.

`make bump VERSION=3.1.0` selects an explicit version. The default patch update
is idempotent against `BASE=origin/main`; repeated invocations do not increment
an already bumped local version. Update the remote tracking reference before
preparing a new push. `BASE=<commit>` supports other already reviewed bases.
The updater validates all fields before writing and restores the original files
if a write fails. It uses the UTC release date to match hosted validation.

Dependency updates, including Dependabot changes to GitHub Actions, require
the same `make bump BASE=origin/main` step before pushing the reviewed update.
The version policy is unchanged: the CI preflight rejects an unchanged version
and explains the command to repair it. No privileged workflow executes code
from an untrusted pull request to modify its metadata automatically.

Benchmark CSV rows include the package version, Git commit, dirty-worktree flag
and actual loaded package path. The temporary installation is removed at exit.
The low-level benchmark script refuses to run without that installation step,
so a previously installed package cannot silently become the benchmark target.
The ensemble benchmark requires every expected extra path to succeed. Runtime
budgets are regression guards; compare versions on the same machine and record
peak memory separately when evaluating allocation changes.

CI distinguishes changes to DESCRIPTION dependency/build fields from a simple
Version/Date update. Package changes retain all OS, minimum-R/testthat,
hard-dependency and coverage checks. Metadata-only and unrelated action updates
still run the Linux package check. Scheduled/manual checks use the full matrix,
and release validation requires successful full checks at the released commit,
dispatching missing checks rather than accepting skipped jobs.
