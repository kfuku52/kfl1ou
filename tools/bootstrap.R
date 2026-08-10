repos <- getOption("repos")
if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@")) {
    repos["CRAN"] <- "https://cloud.r-project.org"
}
options(repos = repos)

if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
}

pak::pkg_install(
    c(
        "deps::.",
        "any::covr",
        "any::rcmdcheck",
        "any::roxygen2"
    ),
    dependencies = TRUE
)
