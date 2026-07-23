# kfl1ou 2.2.1

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
