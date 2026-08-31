# Execute the actual documented R blocks, reducing only the simulation count.
# The ordinary vignette build deliberately skips those expensive blocks.
read_vignette_chunk <- function(path, label) {
    lines <- readLines(path, warn = FALSE)
    start <- grep(paste0("^```\\{r ", label, "([,}])"), lines)
    if (length(start) != 1L || start == length(lines)) {
        stop("Expected one named chunk '", label, "' in ", path, call. = FALSE)
    }
    following <- seq.int(start + 1L, length(lines))
    end <- which(grepl("^```[[:space:]]*$", lines[following]))
    if (!length(end) || end[[1L]] <= 1L) {
        stop("Missing or empty chunk '", label, "' in ", path, call. = FALSE)
    }
    lines[following[seq_len(end[[1L]] - 1L)]]
}

run_documentation_examples <- function(use_installed = FALSE, build_site = FALSE) {
    # A real installation also makes the current package available to PSOCK
    # workers. CI may reuse the checkout just installed by setup-r-dependencies.
    if (!use_installed) {
        if (isNamespaceLoaded("kfl1ou")) {
            stop("Run documentation checks in a fresh Rscript process.", call. = FALSE)
        }
        library <- tempfile("kfl1ou-doc-library-")
        dir.create(library)
        old_paths <- .libPaths()
        on.exit({
            if ("package:kfl1ou" %in% search()) {
                detach("package:kfl1ou", character.only = TRUE, unload = TRUE)
            } else if (isNamespaceLoaded("kfl1ou")) {
                unloadNamespace("kfl1ou")
            }
            .libPaths(old_paths)
            unlink(library, recursive = TRUE)
        }, add = TRUE)
        status <- system2(file.path(R.home("bin"), "R"), c(
            "CMD", "INSTALL", "--preclean", "--no-multiarch", "--with-keep.source",
            paste0("--library=", shQuote(library)), shQuote(normalizePath("."))
        ))
        if (status != 0L) stop("Installing the checkout failed.", call. = FALSE)
        .libPaths(c(library, old_paths))
    }
    loadNamespace("kfl1ou")
    expected <- read.dcf("DESCRIPTION", fields = "Version")[[1L]]
    stopifnot(as.character(getNamespaceVersion("kfl1ou")) == expected)
    if (!use_installed) {
        stopifnot(normalizePath(find.package("kfl1ou")) ==
                  normalizePath(file.path(library, "kfl1ou")))
    }
    message("Documentation examples use kfl1ou ", expected,
            " from ", find.package("kfl1ou"), ".")
    output <- tempfile("kfl1ou-doc-plots-", fileext = ".pdf")
    grDevices::pdf(output)
    on.exit({ grDevices::dev.off(); unlink(output) }, add = TRUE)

    readme <- readLines("README.md", warn = FALSE)
    begin <- which(readme == "<!-- quickstart:start -->")
    end <- which(readme == "<!-- quickstart:end -->")
    stopifnot(length(begin) == 1L, length(end) == 1L, begin < end)
    code <- readme[seq.int(begin + 1L, end - 1L)]
    stopifnot(code[[1L]] == "```r", tail(code, 1L) == "```")
    quickstart <- new.env(parent = globalenv())
    eval(parse(text = code[-c(1L, length(code))]), envir = quickstart)
    stopifnot(inherits(quickstart$fit, "l1ou"), is.finite(quickstart$fit$score))
    message("README quick start passed.")

    run_chunks <- function(file, labels) {
        environment <- new.env(parent = globalenv())
        environment$n_replicates <- 2L
        for (label in labels) {
            message(basename(file), ": ", label)
            code <- read_vignette_chunk(file, label)
            eval(parse(text = code), envir = environment)
        }
        stopifnot(
            environment$support$attempted == environment$n_replicates,
            environment$support$successful + environment$support$failed ==
                environment$n_replicates,
            identical(attr(environment$intervals, "selection"), "full"),
            attr(environment$intervals, "successful") >= 2L
        )
        environment
    }
    run_chunks("vignettes/getting-started.Rmd", c(
        "prepare-data", "bounded-repair", "shift-search", "observation-error",
        "diagnostics", "bootstrap", "full-intervals"
    ))
    multivariate <- run_chunks("vignettes/multivariate-inference.Rmd", c(
        "prepare-data", "fixed-fit", "diagnostics", "covariance-comparison",
        "correlated-search", "bootstrap", "full-intervals"
    ))
    comparison <- multivariate$comparison
    stopifnot(
        identical(comparison$selection, "full"),
        comparison$search.max.nShifts == multivariate$max_shifts,
        comparison$attempted == multivariate$n_replicates,
        comparison$successful >= 1L,
        identical(comparison$null$shift.configuration,
                  multivariate$null_search$shift.configuration),
        identical(comparison$alternative$shift.configuration,
                  multivariate$null_search$shift.configuration)
    )
    message("Documentation examples passed (two simulations per calibration, not for inference).")
    if (build_site) pkgdown::build_site(new_process = FALSE, install = FALSE)
    invisible(TRUE)
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    if (any(!args %in% c("--installed", "--build-site"))) {
        stop("Usage: Rscript tools/check-doc-examples.R [--installed] [--build-site]",
             call. = FALSE)
    }
    run_documentation_examples("--installed" %in% args, "--build-site" %in% args)
}
