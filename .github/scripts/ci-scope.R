ci_description_fields <- function(description) {
    description <- description[setdiff(names(description), c("Version", "Date"))]
    description[order(names(description))]
}

ci_full_checks <- function(paths, description, base.description,
                           event = "pull_request") {
    if (!event %in% c("push", "pull_request")) return(TRUE)
    if (any(grepl(paste0(
        "^(R/|src/|tests/|man/|data/|inst/|vignettes/|tools/|",
        "NAMESPACE$|configure|cleanup|Makefile$|[.]Rbuildignore$|",
        "[.]github/scripts/|[.]github/workflows/(R-CMD-check|extended-checks|release)[.]yaml$)"
    ), paths))) return(TRUE)
    if ("DESCRIPTION" %in% paths) {
        if (is.null(base.description) ||
            !identical(ci_description_fields(description),
                       ci_description_fields(base.description))) return(TRUE)
    }
    FALSE
}

if (sys.nframe() == 0L) {
    base <- Sys.getenv("BASE_SHA", "")
    event <- Sys.getenv("GITHUB_EVENT_NAME", "workflow_dispatch")
    full <- TRUE
    if (nzchar(base) && !grepl("^0+$", base)) {
        paths <- system2("git", c("diff", "--name-only", shQuote(base), "HEAD"),
                         stdout = TRUE)
        previous <- system2("git", c("show", shQuote(paste0(base, ":DESCRIPTION"))),
                            stdout = TRUE)
        if (is.null(attr(paths, "status")) && is.null(attr(previous, "status"))) {
            connection <- textConnection(previous)
            base.description <- read.dcf(connection)[1L, ]
            close(connection)
            full <- ci_full_checks(paths, read.dcf("DESCRIPTION")[1L, ],
                                    base.description, event)
        }
    }
    configurations <- if (full) {
        '[{"os":"macos-latest","r":"release"},{"os":"windows-latest","r":"release"},{"os":"ubuntu-latest","r":"release"}]'
    } else '[{"os":"ubuntu-latest","r":"release"}]'
    output <- c(paste0("full=", tolower(full)), paste0("matrix=", configurations))
    cat(paste(output, collapse = "\n"), "\n")
    destination <- Sys.getenv("GITHUB_OUTPUT", "")
    if (nzchar(destination)) write(output, destination, append = TRUE)
}
