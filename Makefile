R ?= R
RSCRIPT ?= Rscript

.PHONY: bootstrap test document generated check check-full benchmark metadata

bootstrap:
	$(RSCRIPT) tools/bootstrap.R

test:
	$(RSCRIPT) -e 'testthat::test_local(".")'

document:
	$(RSCRIPT) -e 'Rcpp::compileAttributes("."); roxygen2::roxygenise(".")'

generated: document
	git diff --exit-code -- NAMESPACE R/RcppExports.R src/RcppExports.cpp man

check:
	$(RSCRIPT) tools/dev-check.R quick

check-full:
	$(RSCRIPT) tools/dev-check.R full

benchmark:
	$(RSCRIPT) tools/benchmark-search.R

metadata:
	$(RSCRIPT) .github/scripts/check-metadata.R
	$(RSCRIPT) .github/scripts/test-maintenance-scripts.R
