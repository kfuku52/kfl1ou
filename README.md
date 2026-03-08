
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

Now in the shell, after cloning your fork and replacing the asterisks with the correct version number:
```shell
git clone https://github.com/kfuku52/kfl1ou.git
R CMD build kfl1ou
R -e 'install.packages("kfl1ou_*.**.tar.gz")'
```

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
