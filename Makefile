R ?= R
RSCRIPT ?= Rscript
FILTER ?=
VERSION ?= patch
BASE ?= origin/main

.PHONY: bootstrap test document generated doc-examples docs check check-full benchmark metadata bump

bootstrap:
	$(RSCRIPT) tools/bootstrap.R

test:
	$(RSCRIPT) -e 'args <- commandArgs(TRUE); testthat::test_local(".", filter=if (nzchar(args[1])) args[1] else NULL)' '$(FILTER)'

document:
	$(RSCRIPT) -e 'Rcpp::compileAttributes("."); roxygen2::roxygenise(".")'

generated: document
	$(RSCRIPT) tools/check-generated.R

doc-examples:
	$(RSCRIPT) tools/check-doc-examples.R

docs:
	$(RSCRIPT) tools/check-doc-examples.R --build-site
	python3 tools/check_site_links.py docs

check:
	$(RSCRIPT) tools/dev-check.R quick

check-full:
	$(RSCRIPT) tools/dev-check.R full

benchmark:
	$(RSCRIPT) tools/run-benchmark.R

bump:
	BASE_SHA='$(BASE)' $(RSCRIPT) tools/update-metadata.R '$(VERSION)'
	$(RSCRIPT) .github/scripts/check-metadata.R

metadata:
	$(RSCRIPT) .github/scripts/check-metadata.R
	$(RSCRIPT) .github/scripts/test-maintenance-scripts.R
	python3 .github/scripts/test-release-checks.py
	python3 tools/test_doc_tools.py
