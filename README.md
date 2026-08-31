# kfl1ou

[![R-CMD-check](https://github.com/kfuku52/kfl1ou/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/kfuku52/kfl1ou/actions/workflows/R-CMD-check.yaml)
[![CodeQL](https://github.com/kfuku52/kfl1ou/actions/workflows/codeql.yaml/badge.svg?branch=main)](https://github.com/kfuku52/kfl1ou/actions/workflows/codeql.yaml)
[![License: GPL v3+](https://img.shields.io/badge/license-GPL%20v3%2B-blue.svg)](https://github.com/kfuku52/kfl1ou/blob/main/LICENSE)

`kfl1ou` is an independently maintained fork of `l1ou` for detecting changes in
adaptive trait optima on phylogenetic trees with Ornstein–Uhlenbeck models.
It supports convergent regimes, multivariate evolutionary covariance,
observation error, and uncertainty in shifts and phylogenies.

## Install

Requires **R >= 4.1.0** and tools for compiling C/C++ source packages:
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, Xcode
Command Line Tools on macOS, or R development tools and a compiler on Linux.

```r
install.packages("pak")
pak::pkg_install("kfuku52/kfl1ou")
```

Without pak, install the runtime dependencies in R:

```r
install.packages(c("ape", "Rcpp"))
```

Then run in a shell:

```sh
git clone https://github.com/kfuku52/kfl1ou.git
R CMD build --no-build-vignettes kfl1ou
R CMD INSTALL kfl1ou_*.tar.gz
```

This route omits local vignettes; the [guides are available online](https://kfuku52.github.io/kfl1ou/articles/index.html).
To include local vignettes, install `knitr` and `rmarkdown` in R and
[Pandoc](https://pandoc.org/installing.html), then omit `--no-build-vignettes`.

## Quick start

<!-- quickstart:start -->
```r
library(kfl1ou)
data(lizard.tree, lizard.traits)

keep <- lizard.tree$tip.label[seq_len(12)]
tree <- ape::drop.tip(lizard.tree, setdiff(lizard.tree$tip.label, keep))
dat <- adjust_data(tree, lizard.traits[keep, 1, drop = FALSE], quietly = TRUE)
fit <- estimate_shift_configuration(
  dat$tree, dat$Y, criterion = "BIC", max.nShifts = 2
)
plot(fit)
fit$search.diagnostics
```
<!-- quickstart:end -->

`adjust_data()` **normalizes branch lengths to tree height 1 by default**;
use `normalize = FALSE` to retain the original time unit. Tree repair is
separately opt-in. Fully missing tips and invariant traits are removed by
default. See the [data preparation guide](https://kfuku52.github.io/kfl1ou/articles/getting-started.html)
for these controls and conversion back to original time units.

Shift locations describe phylogenetically consistent trait clusters; extant
traits generally do not uniquely identify historical shift times or positions.
Check model adequacy, parameter bounds, search coverage, and sensitivity before
interpreting a selected model.

## Documentation and development

- [Getting started: data, observation error, search and uncertainty](https://kfuku52.github.io/kfl1ou/articles/getting-started.html)
- [Correlated multivariate inference](https://kfuku52.github.io/kfl1ou/articles/multivariate-inference.html)
- [Function reference](https://kfuku52.github.io/kfl1ou/reference/index.html) and [changelog](https://kfuku52.github.io/kfl1ou/news/index.html)
- [Contributing](https://github.com/kfuku52/kfl1ou/blob/main/CONTRIBUTING.md) and [development commands](https://github.com/kfuku52/kfl1ou/blob/main/tools/README.md)

For a development checkout, run `make bootstrap` once, then use focused tests
and `make check`; run `make check-full` before submission. These guides and the
reference site are the maintained documentation; a separate Wiki is not used.

## Citation and license

Run `citation("kfl1ou")` for the package citation and the original l1ou article;
please cite both. The fork is distributed under **GPL-3.0-or-later**.
See [LICENSE](https://github.com/kfuku52/kfl1ou/blob/main/LICENSE) for the license
and [NOTICE](https://github.com/kfuku52/kfl1ou/blob/main/NOTICE) for provenance.
