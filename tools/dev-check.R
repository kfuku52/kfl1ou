args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "quick"
if (!mode %in% c("quick", "full")) {
    stop("Usage: Rscript tools/dev-check.R [quick|full]", call. = FALSE)
}
if (!requireNamespace("rcmdcheck", quietly = TRUE)) {
    stop("Install development dependencies first with `make bootstrap`.",
         call. = FALSE)
}

repository <- normalizePath(".", mustWork = TRUE)
build_dir <- tempfile("kfl1ou-check-")
dir.create(build_dir)
on.exit(unlink(build_dir, recursive = TRUE, force = TRUE), add = TRUE)

old_dir <- setwd(build_dir)
on.exit(setwd(old_dir), add = TRUE)
build_status <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "build", shQuote(repository), "--compact-vignettes=gs+qpdf")
)
if (!identical(build_status, 0L)) {
    stop("R CMD build failed.", call. = FALSE)
}

tarballs <- list.files(build_dir, pattern = "^kfl1ou_.*[.]tar[.]gz$", full.names = TRUE)
if (length(tarballs) != 1L) {
    stop("Expected exactly one source package after R CMD build.", call. = FALSE)
}

check_args <- if (identical(mode, "full")) "--as-cran" else "--no-manual"
rcmdcheck::rcmdcheck(tarballs[[1L]], args = check_args, error_on = "warning")
