metadata_value <- function(path, key) {
    lines <- readLines(path, warn = FALSE)
    match <- grep(paste0("^", key, ":[[:space:]]*"), lines, value = TRUE)
    if (length(match) != 1L) {
        stop(sprintf("Expected one `%s` field in %s.", key, path), call. = FALSE)
    }
    value <- sub(paste0("^", key, ":[[:space:]]*"), "", match)
    sub('^(["\' ]?)(.*)\\1$', "\\2", trimws(value))
}

metadata_fail <- function(message) {
    stop(paste("Metadata consistency check failed:", message), call. = FALSE)
}

check_metadata <- function(root = ".", base_sha = Sys.getenv("BASE_SHA", "")) {
    description_path <- file.path(root, "DESCRIPTION")
    description <- read.dcf(description_path)[1L, ]
    version <- unname(description[["Version"]])
    package_date <- as.Date(unname(description[["Date"]]))
    if (is.na(package_date) || package_date > Sys.Date()) {
        metadata_fail("DESCRIPTION Date must be valid and must not be in the future.")
    }

    citation_version <- metadata_value(file.path(root, "CITATION.cff"), "version")
    if (!identical(version, citation_version)) {
        metadata_fail(sprintf("DESCRIPTION is %s but CITATION.cff is %s.",
                              version, citation_version))
    }

    news <- readLines(file.path(root, "NEWS.md"), warn = FALSE)
    first_heading <- grep("^# kfl1ou ", news, value = TRUE)[1L]
    expected_heading <- paste("# kfl1ou", version)
    if (is.na(first_heading) || !identical(first_heading, expected_heading)) {
        metadata_fail(sprintf("the first NEWS heading must be `%s`.", expected_heading))
    }

    major <- strsplit(version, ".", fixed = TRUE)[[1L]][1L]
    security <- paste(readLines(file.path(root, "SECURITY.md"), warn = FALSE),
                      collapse = "\n")
    if (!grepl(paste0("current ", major, "[.]x release line"), security)) {
        metadata_fail(sprintf("SECURITY.md must identify the current %s.x release line.",
                              major))
    }

    base_sha <- trimws(base_sha)
    if (nzchar(base_sha) && !grepl("^0+$", base_sha)) {
        git_root <- normalizePath(root, mustWork = TRUE)
        valid_sha <- system2("git", c("-C", shQuote(git_root), "cat-file", "-e",
                                       shQuote(paste0(base_sha, "^{commit}"))),
                             stdout = FALSE, stderr = FALSE)
        if (identical(valid_sha, 0L)) {
            old_description <- tryCatch(
                system2("git", c("-C", shQuote(git_root), "show",
                                  shQuote(paste0(base_sha, ":DESCRIPTION"))),
                        stdout = TRUE, stderr = FALSE),
                error = function(error) character()
            )
            if (length(old_description)) {
                old_file <- tempfile("description-")
                on.exit(unlink(old_file), add = TRUE)
                writeLines(old_description, old_file)
                old_version <- unname(read.dcf(old_file)[1L, "Version"])
                if (identical(version, old_version)) {
                    metadata_fail(sprintf(
                        "Version %s is unchanged from BASE_SHA %s.", version, base_sha
                    ))
                }
            }
        }
    }

    message(sprintf("Metadata is consistent for kfl1ou %s (%s).", version, package_date))
    invisible(TRUE)
}

if (sys.nframe() == 0L) {
    root <- Sys.getenv("KFL1OU_METADATA_ROOT", ".")
    check_metadata(root)
}
