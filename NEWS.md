# kfl1ou 2.5.4

- Keep CodeQL initialization and analysis on the same pinned release, and
  group their future Dependabot updates into one pull request.

# kfl1ou 2.5.3

- Accelerate full multivariate OU fits with direct covariance construction,
  cached tree geometry, cost- and memory-aware automatic likelihood selection,
  vectorized pruning preprocessing and diagonal message algebra, and optional
  parallel execution of deterministic optimizer starts within the configured
  CPU budget.
- Extend the deterministic benchmark gate to cover automatic full-covariance
  fits and tighten the existing exhaustive and ensemble search budgets.

# kfl1ou 2.5.1

- Keep the CodeQL initialization and analysis actions on the same pinned
  release so dependency updates cannot create a mixed-version analysis job.

# kfl1ou 2.5.0

- Correct convergent-regime design matrices for nested shifts, refit the final
  design at the estimated adaptation rate, retain continuous OU weights when
  alpha is fixed, and base BIC/AICc penalties on the dimensions of the models
  actually fitted.
- Preserve finite expected mean shifts at the exact Brownian boundary while
  marking unidentifiable optimum shifts as missing, propagate fixed-alpha
  bounds consistently, and handle observation-error entries for unobserved
  multivariate cells during latent simulation.
- Strengthen tree, shift, candidate-edge, replicate, ensemble-weight, and
  simulation validation; reconstruct complete latent fitted means for
  missing-data simulations and distinguish raw from trait-whitened diagnostic
  correlations.
- Label non-nested and regularized trait-covariance comparisons explicitly,
  calibrate them only when justified, and report their actual fitted parameter
  difference.
- Harden native pruning, covariance, configuration-score, and group-lasso
  routines against interrupted unwinding, size overflow, input mutation,
  subnormal branch loss, non-finite scores, and unstable large-scale arithmetic.
- Expand public-API and regression coverage, execute a lightweight
  multivariate vignette, test the minimum supported testthat and reduced
  dependency sets, build the PDF manual, audit the R dependency closure with
  CycloneDX and OSV, and gate releases on successful checks and provenance
  attestation.

# kfl1ou 2.4.3

- Stabilize OU branch transformations and design matrices near the Brownian
  motion boundary, and warn when a stationary-root covariance is intrinsically
  ill-conditioned at extremely small adaptation rates.
- Ignore machine-precision group-lasso coefficients when identifying active
  supports and add randomized objective-equivalence regression coverage for
  the C++ and vendored solvers.
- Run deterministic search benchmarks for pushes and pull requests, enforce
  generous performance-regression budgets, and use stable `log1p` evaluations
  in native and R likelihood calculations.
- Derive the installed package citation version from `DESCRIPTION` so release
  metadata cannot silently drift.

# kfl1ou 2.4.2

- Validate native tree topology, dimensions, branch lengths, group-lasso
  controls, and numerical inputs before entering allocation and traversal
  kernels; use R-managed memory so native errors cannot leak workspaces.
- Stop regularized convergent-regime searches once all shifts have merged into
  one regime instead of attempting to build an empty relation matrix.
- Replace unsafe sequence construction in legacy loops and add regression
  coverage for malformed native inputs, multi-shift convergence, and all plot
  annotation paths.

# kfl1ou 2.4.1

- Prevent integer overflow when sizing native three-point likelihood
  workspaces, reject infeasible allocations before entering C, and handle
  allocation failure without dereferencing null pointers.

# kfl1ou 2.4.0

- Exercise the requested sparse-search backends explicitly in regression tests
  instead of allowing automatic exhaustive search to bypass them.
- Harden continuous integration with pinned actions, concurrency cancellation,
  timeouts, generated-file drift checks, patch and critical-file coverage
  gates, dependency review, CodeQL, and scheduled extended platform checks.
- Add scheduled AddressSanitizer, UndefinedBehaviorSanitizer, and Valgrind
  checks using the CRAN-like R-hub containers.
- Add reproducible performance and determinism benchmarks for exact and
  ensemble shift searches.
- Add Dependabot configuration, release automation, a pkgdown site workflow,
  contribution and security policies, and machine-readable citation metadata.
- Publish coverage and benchmark artifacts so regressions can be inspected
  without rerunning a workflow.

# kfl1ou 2.3.0

- Add automatically certified exhaustive shift searches for small
  configuration spaces, seeded ensemble-path augmentation for larger searches,
  candidate-space coverage diagnostics, and complete candidate profiles.
- Make full-search bootstrap intervals work for models created by `fit_OU()`;
  non-selected exact-edge coefficients now contribute zero, with exact-edge
  and tip-partition reselection probabilities reported explicitly.
- Add canonical tip partitions, equivalent-scenario enumeration, partition
  bootstrap summaries, weighted tip co-assignment, and tree-ensemble
  sensitivity analysis.
- Add projected residual contrasts, fitted-model predictive checks,
  conditional tip outlier scores, clade diffusion-rate diagnostics, and a
  scan-calibrated clade rate-shift likelihood comparison.
- Retain trait-specific branch effects in multivariate sparse paths by default,
  iteratively regenerate full-covariance candidates after mean refits, and add
  a localization-aware multivariate pBIC sensitivity criterion.
- Replace the fixed covariance-shrinkage heuristic with an analytic estimate
  based on phylogenetically whitened residuals, while retaining an explicit
  fixed-lambda option.
- Add replicated-measurement fitting, measurement-error profiles, full-pipeline
  covariance bootstrap selection, and root/alpha/criterion sensitivity grids.
- Restore Windows builds by isolating POSIX-only runtime symbol lookup, and
  harden the native three-point likelihood boundary against malformed input.
- Keep non-evaluated vignette examples out of extracted check code, expand
  residual-bootstrap regression coverage, and enforce coverage floors in CI.
- Prevent Wald intervals from interpreting covariance parameters as fixed
  adaptation rates; penalized Hessians are no longer exposed as ordinary
  parameter covariance matrices.
- Correct information-criterion parameter counts when alpha is fixed, and
  support the exact alpha-zero Brownian boundary in fixed-alpha fits and
  likelihood profiles.
- Restore warning and random-number state on all error and seeded simulation
  paths, and make parametric full-search shift bootstrap the inferential
  default with explicit failure accounting.
- Correct the native root-state traversal and safely support root polytomies.
- Use a single shrinkage-penalty definition across dense, pruning, complete,
  and incomplete fits, and honor explicit likelihood-engine requests.
- Remove dense trait-by-tip covariance construction from pruning
  initialization and add direct Gaussian tree simulation.
- Validate numeric traits, finite values, unique labels and positive branch
  lengths at public boundaries, including strict integer/logical controls and
  known-error dimensions; retain normalized-tree time scales and provide
  original-unit and half-life conversion helpers.
- Separate coefficient `vcov()` semantics from evolutionary innovation
  covariance through `evolutionary_vcov()`.
- Add alpha likelihood profiles, clade-level bootstrap summaries, Monte Carlo
  error for covariance comparisons, and safer refitted model averaging.

- Generalize the correlated multivariate OU model to a diagonal drift matrix,
  allowing a separate adaptation rate for every trait together with a full
  evolutionary innovation covariance matrix.
- Add an exact Gaussian tree-pruning likelihood for diagonal-drift OU models.
  It supports fixed or stationary roots, missing trait entries, known
  observation variances, and estimated measurement error; dense and pruning
  likelihoods are regression-tested against one another.
- Add deterministic multi-start optimization, boundary checks, covariance
  conditioning diagnostics, optional numerical Hessians, and recorded
  optimizer runs.
- Add covariance shrinkage toward a diagonal target for high-dimensional or
  poorly conditioned trait sets.
- Add standard model methods (`coef`, `fitted`, `residuals`, `logLik`, `nobs`,
  `vcov`, `predict`, `simulate`, and `confint`) plus `diagnose_l1ou()`.
- Add parametric-bootstrap comparison of diagonal versus full trait
  covariance, bootstrap shift co-selection summaries, and
  information-criterion model averaging across shift configurations.

# kfl1ou 2.1.0

- Add an opt-in, joint multivariate OU likelihood with a shared adaptation rate
  and an estimated positive-definite evolutionary trait covariance matrix.
  The full model supports trait-specific missing data, observation error,
  convergent-regime refitting, BIC model selection, and joint bootstrap
  simulation.
- Keep the historical block-diagonal, trait-specific-alpha likelihood as
  `trait.covariance = "diagonal"` for backward compatibility.
- Reject information-criterion penalties that have not been derived for the
  correlated multivariate model instead of applying them heuristically.

# kfl1ou 2.0.3

- Correct the likelihood for known tip-specific observation variances so that
  phylogenetic and observation variances are scaled separately.
- Refit convergent-regime models under their equality constraints and return
  constrained coefficients, fitted values, residuals, optima, variances, and
  likelihoods. Convergent fits now preserve the selected root model.
- Use the observed marginal covariance for multivariate GLS whitening and make
  the multivariate residual bootstrap support trait-specific missingness.
- Continue scanning non-monotone sparse solution paths after configurations
  that exceed `max.nShifts`.
- Make `alpha_upper_bound()` invariant to tip numbering and restore its
  documented minimum-external-branch definition.
- Make automatic tree repair opt-in and enforce the requested relative-change
  bound for every repair backend.
- Add regression tests, declared optional dependencies, CI, and explicit
  documentation of the statistical limitations of OU shift inference.
