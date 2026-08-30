args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "quick"
if (!mode %in% c("quick", "full")) {
    stop("Usage: Rscript tools/dev-check.R [quick|full]", call. = FALSE)
}
if (!requireNamespace("rcmdcheck", quietly = TRUE)) {
    stop("Install development dependencies first with `make bootstrap`.",
         call. = FALSE)
}

run_check <- function(mode) {
    repository <- normalizePath(".", mustWork = TRUE)
    output <- Sys.getenv("KFL1OU_CHECK_OUTPUT", "")
    if (nzchar(output)) {
        dir.create(output, recursive = TRUE, showWarnings = FALSE)
        output <- normalizePath(output, mustWork = TRUE)
    }
    build_dir <- tempfile("kfl1ou-check-", tmpdir = if (nzchar(output)) output else tempdir())
    dir.create(build_dir)
    old_dir <- setwd(build_dir)
    on.exit({
        setwd(old_dir)
        if (!nzchar(output)) unlink(build_dir, recursive = TRUE, force = TRUE)
    }, add = TRUE)
    build_args <- if (identical(mode, "full")) "--compact-vignettes=gs+qpdf" else
        c("--no-build-vignettes", "--no-manual")
    build_status <- system2(
        file.path(R.home("bin"), "R"),
        c("CMD", "build", shQuote(repository), build_args)
    )
    if (!identical(build_status, 0L)) {
        stop("R CMD build failed.", call. = FALSE)
    }

    tarballs <- list.files(build_dir, pattern = "^kfl1ou_.*[.]tar[.]gz$", full.names = TRUE)
    if (length(tarballs) != 1L) {
        stop("Expected exactly one source package after R CMD build.", call. = FALSE)
    }

    check_args <- if (identical(mode, "full")) "--as-cran" else
        c("--no-manual", "--ignore-vignettes")
    rcmdcheck::rcmdcheck(tarballs[[1L]], args = check_args, error_on = "warning",
                         check_dir = file.path(build_dir, "check"))
}

run_check(mode)
