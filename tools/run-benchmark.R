run_benchmark <- function() {
    root <- normalizePath(".", mustWork = TRUE)
    library <- tempfile("kfl1ou-benchmark-library-")
    dir.create(library)
    on.exit(unlink(library, recursive = TRUE), add = TRUE)
    status <- system2(file.path(R.home("bin"), "R"), c(
        "CMD", "INSTALL", "--preclean", "--no-multiarch", "--with-keep.source",
        paste0("--library=", shQuote(library)), shQuote(root)
    ))
    if (!identical(status, 0L)) stop("Installing the checkout failed.", call. = FALSE)
    old.paths <- .libPaths()
    on.exit(.libPaths(old.paths), add = TRUE)
    .libPaths(c(library, old.paths))
    loadNamespace("kfl1ou")
    loaded <- normalizePath(find.package("kfl1ou"))
    stopifnot(identical(loaded, normalizePath(file.path(library, "kfl1ou"))))
    sha <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
    dirty <- length(system2("git", c("status", "--porcelain"), stdout = TRUE)) > 0L
    options(kfl1ou.benchmark.source = list(
        package.version = as.character(utils::packageVersion("kfl1ou")),
        git.sha = sha, worktree.dirty = dirty, loaded.from = loaded
    ))
    source(file.path(root, "tools", "benchmark-search.R"), local = new.env())
}

run_benchmark()
