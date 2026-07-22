
### Tools for detecting past changes in the expected mean trait values and studying trait evolution from comparative data

The kfl1ou package is an independently maintained fork of l1ou. It provides functions to study trait evolution from comparative data and detect past changes in the expected mean trait values, as well as convergent evolution. It uses the Ornstein-Uhlenbeck process along a phylogenetic tree, which can model a changing adaptive landscape over time and over lineages. 

### Origin and license

`kfl1ou` is an independently maintained derivative of `l1ou`.
This package is distributed under the GNU General Public License, version 3 or
later (`GPL (>= 3)` in `DESCRIPTION`). The full license text is available in
[LICENSE](LICENSE), and provenance notes for this fork are summarized in
[NOTICE](NOTICE).

### Install using the devtools package
From within R:
```r
install.packages("devtools")
library(devtools)
install_github("kfuku52/kfl1ou")
```
Windows users will first need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

### Install without the devtools package

To resolve hard dependencies, first install the following packages from CRAN.
From within R:
```r
install.packages(c("ape", "Rcpp"))
```

Now in the shell, after cloning the repository:
```shell
git clone https://github.com/kfuku52/kfl1ou.git
R CMD build kfl1ou
R CMD INSTALL kfl1ou_*.tar.gz
```

### Data preparation and tree repair

`adjust_data()` reorders and validates data but does not alter branch lengths by
default. If a nearly ultrametric tree needs repair, request it explicitly and
set the maximum acceptable relative branch-length change:

```r
prepared <- adjust_data(
  tree,
  Y,
  repair.tree = TRUE,
  ultrametric.tolerance = 0.01
)
```

Repair fails instead of silently exceeding the supplied tolerance. Review any
reported branch-length change before interpreting fitted shift locations.

### Observation error handling

`kfl1ou` supports two observation-error inputs:

- `input_error`: known tip-specific observation error variances.
- `measurement_error = TRUE`: estimate an additional shared observation-error variance `sigma2_error` for each trait.

These can be used separately or together. When both are supplied, the fitted covariance is:

`Sigma_phylo + Diag(input_error) + sigma2_error I`

Typical usage from within R:

```r
## known tip-specific observation error only
fit_known <- estimate_shift_configuration(
  tree,
  Y,
  input_error = known_tip_variances
)

## estimate an additional shared measurement-error variance
fit_estimated <- estimate_shift_configuration(
  tree,
  Y,
  measurement_error = TRUE
)

## use known tip-specific error and estimate extra unexplained error
fit_joint <- estimate_shift_configuration(
  tree,
  Y,
  input_error = known_tip_variances,
  measurement_error = TRUE
)
```

For multivariate data, `input_error` can be a matrix with one column per trait.

### Correlated multivariate OU model

Set `trait.covariance = "full"` to estimate evolutionary covariance among
traits instead of fitting a block-diagonal trait likelihood:

```r
fit_correlated <- estimate_shift_configuration(
  tree,
  Y,
  max.nShifts = 2,
  criterion = "BIC",
  trait.covariance = "full"
)

fit_correlated$trait.covariance
fit_correlated$trait.correlation
fit_correlated$tip.trait.correlation
fit_correlated$joint.logLik
evolutionary_vcov(fit_correlated)
```

By default this mode fits the separable covariance `Omega %x% C(alpha)`, where
`Omega` is a positive-definite evolutionary innovation covariance matrix and
`C(alpha)` is the OU covariance among tips. To estimate a distinct adaptation
rate for every trait while retaining correlated innovations, use:

```r
fit_general <- fit_OU(
  tree,
  Y,
  shift.configuration = integer(),
  criterion = "BIC",
  trait.covariance = "full",
  alpha.structure = "diagonal",
  likelihood.engine = "auto"
)
```

The diagonal-alpha model uses a full innovation covariance and the exact
cross-trait covariance induced by the trait-specific OU rates. Both dense and
Gaussian tree-pruning likelihood engines support trait-specific missing
entries, known observation variances, estimated measurement-error variances,
and fixed or stationary roots. `likelihood.engine = "auto"` retains the
matrix-normal profile for complete, unregularized shared-alpha data and selects
pruning for larger non-separable problems. Explicit `"dense"` and `"pruning"`
requests are honored. Pruning initialization and tree-based simulation do not
construct the full trait-by-tip covariance matrix.

`trait.covariance` and `trait.correlation` describe instantaneous evolutionary
innovations. When alpha differs among traits, `tip.trait.covariance` and
`tip.trait.correlation` report the induced marginal residual association at an
ultrametric tip; the distinction is important because the two correlations are
then not generally equal.

For many traits or few residual contrasts, use
`covariance.regularization = "shrinkage"`; an explicit
`regularization.lambda` can be supplied or selected automatically. Full
covariance requires `criterion = "BIC"`: the package deliberately
rejects pBIC, pBICess, mBIC, and AICc because their penalties have not been
derived for this joint model. The historical `trait.covariance = "diagonal"`
mode remains the default and retains trait-specific alpha estimates. Shrinkage
has the same correlation-penalty interpretation for dense and pruning engines,
including incomplete data. BIC evaluated at a regularized estimate is reported
as a sensitivity score rather than a calibrated marginal-likelihood
approximation.

Numerical and inferential checks are available directly from fitted objects:

```r
diagnose_l1ou(fit_general)
confint(
  fit_general, method = "parametric", selection = "full",
  nsim = 200, seed = 1
)
compare_trait_covariance(tree, Y, nboot = 200, seed = 2)

support <- l1ou_bootstrap_support(
  fit_correlated, nItrs = 200, type = "parametric", seed = 3
)
summarize_shift_uncertainty(support, tree)
averaged <- model_average_l1ou(fit_correlated, delta.max = 10)
```

### Interpretation and statistical limitations

- Shift configurations describe phylogenetically consistent clusters. Tip-only
  data generally cannot identify the exact number, branch position, or timing of
  historical shifts without additional information such as fossils.
- In the default diagonal covariance mode, multivariate fits treat trait
  residuals as independent. Full covariance with `alpha.structure =
  "diagonal"` handles correlated innovations and trait-specific OU rates. A
  fully non-diagonal drift matrix (directional coupling among traits) is still
  outside the implemented model.
- OU optima and adaptation rates can be weakly identified even with many tips.
  Treat estimates at parameter bounds as sensitivity warnings and compare
  biologically plausible bounds and root models. `profile_alpha_l1ou()` and the
  `alpha_to_half_life()` helpers make these checks easier to report.
- pBIC for unconstrained shifts follows the derivation in Khabbazian et al.
  (2016). Its use for convergent equality constraints remains heuristic; report
  sensitivity to AICc/BIC and bootstrap support rather than treating one
  criterion as calibrated certainty.
- `adjust_data(normalize = TRUE)` stores the original tree-height factor.
  `original_time_parameters()` converts alpha and evolutionary variance rates
  back to the original branch-length time unit.
- `vcov()` follows the R convention of coefficient sampling covariance and may
  return an explicitly unavailable matrix. Use `evolutionary_vcov()` for the
  fitted evolutionary innovation covariance.

Convergent-regime fits are refitted under their equality constraints. Returned
coefficients, fitted values, residuals, optima, likelihoods, and variance
parameters therefore all correspond to the constrained model.
