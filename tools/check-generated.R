check_generated <- function(root = ".") {
    paths <- c("NAMESPACE", "R/RcppExports.R", "src/RcppExports.cpp", "man")
    status <- system2(
        "git",
        c("-C", shQuote(normalizePath(root, mustWork = TRUE)),
          "status", "--porcelain=v1", "--untracked-files=all", "--", paths),
        stdout = TRUE, stderr = TRUE
    )
    exit_status <- attr(status, "status")
    if (!is.null(exit_status) && exit_status != 0L) {
        stop("Could not inspect generated files: ", paste(status, collapse = "\n"),
             call. = FALSE)
    }
    if (length(status)) {
        stop(
            "Generated files differ from HEAD (including staged or untracked files):\n",
            paste(status, collapse = "\n"),
            "\nRun `make document`, review and commit its output, then run `make generated`.",
            call. = FALSE
        )
    }
    message("Generated files match HEAD.")
    invisible(TRUE)
}

if (sys.nframe() == 0L) check_generated()
