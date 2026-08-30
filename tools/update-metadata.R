replace_metadata_field <- function(lines, key, value) {
    index <- grep(paste0("^", key, ":"), lines)
    if (length(index) != 1L) stop("Expected one ", key, " field.", call. = FALSE)
    lines[index] <- paste0(key, ": ", value)
    lines
}

update_metadata <- function(root = ".", version = "patch", base.version = NULL,
                            note = "Maintenance and compatibility updates.",
                            date = as.Date(Sys.time(), tz = "UTC")) {
    paths <- c("DESCRIPTION", "CITATION.cff", "NEWS.md", "SECURITY.md")
    original <- lapply(file.path(root, paths), readLines, warn = FALSE)
    names(original) <- paths
    current <- unname(read.dcf(file.path(root, "DESCRIPTION"))[1L, "Version"])
    valid <- function(x) length(x) == 1L && !is.na(x) &&
        grepl("^[0-9]+[.][0-9]+[.][0-9]+$", x)
    if (!valid(current) || (!is.null(base.version) && !valid(base.version))) {
        stop("Current and base versions must have the form MAJOR.MINOR.PATCH.",
             call. = FALSE)
    }
    if (identical(version, "patch")) {
        # Idempotent against a base commit: repeated local edits need only one
        # bump per push. This also works for checked-out dependency updates.
        if (!is.null(base.version) && utils::compareVersion(current, base.version) > 0L) {
            version <- current
        } else {
            baseline <- if (is.null(base.version)) current else base.version
            parts <- as.integer(strsplit(baseline, ".", fixed = TRUE)[[1L]])
            parts[3L] <- parts[3L] + 1L
            version <- paste(parts, collapse = ".")
        }
    }
    if (!valid(version) || utils::compareVersion(version, current) < 0L ||
        (!is.null(base.version) && utils::compareVersion(version, base.version) <= 0L)) {
        stop("The target version must not decrease and must exceed the base version.",
             call. = FALSE)
    }
    if (length(date) != 1L || is.na(date) || date > as.Date(Sys.time(), tz = "UTC")) {
        stop("The release date must be a valid date no later than today in UTC.",
             call. = FALSE)
    }
    updated <- original
    updated$DESCRIPTION <- replace_metadata_field(
        replace_metadata_field(updated$DESCRIPTION, "Version", version),
        "Date", format(date)
    )
    updated$CITATION.cff <- replace_metadata_field(updated$CITATION.cff, "version", version)
    heading <- paste("# kfl1ou", version)
    first <- grep("^# kfl1ou ", updated$NEWS.md, value = TRUE)[1L]
    if (!identical(first, heading)) {
        updated$NEWS.md <- c(heading, "", paste0("- ", note), "", updated$NEWS.md)
    }
    major <- strsplit(version, ".", fixed = TRUE)[[1L]][1L]
    if (sum(grepl("current [0-9]+[.]x release line", updated$SECURITY.md)) != 1L) {
        stop("SECURITY.md must contain one current release line.", call. = FALSE)
    }
    updated$SECURITY.md <- sub("current [0-9]+[.]x release line",
                              paste0("current ", major, ".x release line"),
                              updated$SECURITY.md)
    changed <- !vapply(seq_along(paths), function(i) {
        identical(original[[i]], updated[[i]])
    }, logical(1))
    # Validate everything before any write and restore all originals if a
    # write fails, so a partial update cannot leave mixed release metadata.
    tryCatch({
        for (i in which(changed)) writeLines(updated[[i]], file.path(root, paths[i]))
    }, error = function(error) {
        for (i in which(changed)) writeLines(original[[i]], file.path(root, paths[i]))
        stop(error)
    })
    message("Release metadata synchronized: ", version, " (", format(date), ").")
    invisible(version)
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    version <- if (length(args)) args[[1L]] else "patch"
    root <- Sys.getenv("KFL1OU_METADATA_ROOT", ".")
    base <- Sys.getenv("BASE_SHA", "")
    base.version <- NULL
    if (nzchar(base)) {
        lines <- system2("git", c("-C", shQuote(normalizePath(root)), "show",
                                  shQuote(paste0(base, ":DESCRIPTION"))), stdout = TRUE)
        if (!is.null(attr(lines, "status"))) stop("Cannot read BASE_SHA.", call. = FALSE)
        connection <- textConnection(lines)
        base.version <- unname(read.dcf(connection)[1L, "Version"])
        close(connection)
    }
    update_metadata(root, version, base.version)
}
