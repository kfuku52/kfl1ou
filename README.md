
### Tools for detecting past changes in the expected mean trait values and studying trait evolution from comparative data

The kfl1ou package is an independently maintained fork of l1ou. It provides functions to study trait evolution from comparative data and detect past changes in the expected mean trait values, as well as convergent evolution. It uses the Ornstein-Uhlenbeck process along a phylogenetic tree, which can model a changing adaptive landscape over time and over lineages. 
<!--Detection of evolutionary shifts in trait evolution from extant taxa is motivated by the study of convergent evolution, or to correlate shifts in traits with habitat changes or with changes in other phenotypes.-->
Shifts can be detected from multiple traits, assuming that all traits shifted along the same lineages. Estimation is very fast thanks to lasso techniques, and the user can choose from multiple information criteria for model selection, including a phylognetic BIC. 
For backward compatibility, fitted objects still use class `l1ou` and retain the component name `l1ou.options`.
The `nCores` argument is treated as a total CPU budget for kfl1ou: when process parallelism is used, nested kfl1ou calls run sequentially and BLAS/OpenMP threads are reduced per process when supported so the full run respects that budget.
Citation: 

- M. Khabbazian, R. Kriebel, K. Rohe, and C&eacute;cile An&eacute;.
Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models.
Methods in Ecology and Evolution, 7(7):811–824.
[doi:10.1111/2041-210X.12534](http://dx.doi.org/10.1111/2041-210X.12534)

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
install_github("glmgen/genlasso")
install_github("kfuku52/kfl1ou")
```
Windows users will first need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

### Install without the devtools package

To resolve dependencies, first install the following packages from CRAN, then the knitr package.
From within R:
```r
install.packages(c("igraph", "phylolm", "magic", "Rcpp"))
install.packages("knitr")
```

Download genlasso version 1.3 R package from CRAN archive [(link)](https://cran.r-project.org/src/contrib/Archive/genlasso/genlasso_1.3.tar.gz). 
From within R:
```r
 install.packages("genlasso_1.3.tar.gz", repos=NULL, type="source")
 ```

Now in the shell, after cloning your fork and replacing the asterisks with the correct version number:
```shell
git clone https://github.com/kfuku52/kfl1ou.git
R CMD build kfl1ou
R -e 'install.packages("kfl1ou_*.**.tar.gz")'
```

### Version notes 

major changes are indicated below.

- v2.0.0 (2026-03-07): renamed package to kfl1ou and continued independent maintenance
- v1.43 (2022-08-05): compatibility updates
- v1.42 (2019-02-10): bug fix
- v1.41 (2017-06-18): small bug fixes
- v1.40 (2017-03-27):
  * intercept correctly handled after noise-whitening
  (results may change for variables with a mean far from 0)
  * bug fix in the function calculating the square-root (and inverse) of the
  phylogenetic covariance matrix.
- v1.25: the penalty term in the AICc score now considers the intercept as a free variable.
  The change only affects the final value of the AICc score.
- v1.23: "estimate\_convergent\_regimes" function accepts multiple traits. 
- v1.22: 
	* the scores returned by "estimate\_shift\_configuration” function are now for the non-normalized, original data.
	* "estimate\_shift\_configuration" function also accepts multiple traits with missing values. 
