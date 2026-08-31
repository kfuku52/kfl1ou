source(file.path(".github", "scripts", "check-metadata.R"), local = TRUE)
source(file.path(".github", "scripts", "ci-helpers.R"), local = TRUE)
source(file.path(".github", "scripts", "ci-scope.R"), local = TRUE)
source(file.path("tools", "update-metadata.R"), local = TRUE)
source(file.path("tools", "check-generated.R"), local = TRUE)

fixture <- tempfile("kfl1ou-metadata-")
dir.create(fixture)
on.exit(unlink(fixture, recursive = TRUE, force = TRUE), add = TRUE)

writeLines(c(
    "Package: kfl1ou",
    "Version: 9.8.7",
    paste("Date:", as.Date(Sys.time(), tz = "UTC"))
), file.path(fixture, "DESCRIPTION"))
writeLines("version: 9.8.7", file.path(fixture, "CITATION.cff"))
writeLines("# kfl1ou 9.8.7", file.path(fixture, "NEWS.md"))
writeLines("Security fixes target the current 9.x release line.",
           file.path(fixture, "SECURITY.md"))
stopifnot(isTRUE(check_metadata(fixture, base_sha = "")))

writeLines("version: 9.8.6", file.path(fixture, "CITATION.cff"))
bad_citation <- tryCatch({
    check_metadata(fixture, base_sha = "")
    FALSE
}, error = function(error) grepl("CITATION", conditionMessage(error)))
stopifnot(isTRUE(bad_citation))

writeLines("version: 9.8.7", file.path(fixture, "CITATION.cff"))
writeLines("# kfl1ou 9.8.6", file.path(fixture, "NEWS.md"))
bad_news <- tryCatch({
    check_metadata(fixture, base_sha = "")
    FALSE
}, error = function(error) grepl("NEWS", conditionMessage(error)))
stopifnot(isTRUE(bad_news))

stopifnot(identical(update_metadata(fixture, base.version = "9.8.7"), "9.8.8"))
stopifnot(isTRUE(check_metadata(fixture, base_sha = "")))
first.update <- lapply(list.files(fixture, full.names = TRUE), readLines)
stopifnot(identical(update_metadata(fixture, base.version = "9.8.7"), "9.8.8"))
stopifnot(identical(first.update, lapply(list.files(fixture, full.names = TRUE), readLines)))
stopifnot(identical(update_metadata(fixture, "10.0.0"), "10.0.0"))
stopifnot(isTRUE(check_metadata(fixture, base_sha = "")))
before.invalid <- readLines(file.path(fixture, "DESCRIPTION"))
invalid.version <- tryCatch({
    update_metadata(fixture, "1.0.0")
    FALSE
}, error = function(error) TRUE)
stopifnot(invalid.version,
          identical(before.invalid, readLines(file.path(fixture, "DESCRIPTION"))))

base.description <- c(Package = "kfl1ou", Version = "3.0.5", Date = "2026-08-18",
                      Imports = "ape (>= 4.0), Rcpp")
description <- base.description
description[c("Version", "Date")] <- c("3.0.6", "2026-08-30")
action.paths <- c("DESCRIPTION", "NEWS.md", "CITATION.cff",
                  ".github/workflows/codeql.yaml")
stopifnot(!ci_full_checks(action.paths, description, base.description))
stopifnot(ci_full_checks(action.paths, description, base.description, "workflow_dispatch"))
stopifnot(ci_full_checks(character(), description, base.description, "schedule"))
stopifnot(ci_full_checks("R/model_api.R", description, base.description))
stopifnot(ci_full_checks(".github/workflows/R-CMD-check.yaml", description, base.description))
description["Imports"] <- "ape (>= 4.1), Rcpp"
stopifnot(ci_full_checks(action.paths, description, base.description))

stopifnot(identical(
    parse.dependencies("R (>= 4.1), ape (>= 4.0), Rcpp, ape"),
    c("ape", "Rcpp")
))
stopifnot(identical(
    description.dependencies(list(Depends="R, ape", Imports="Rcpp")),
    c("ape", "Rcpp")
))
stopifnot(identical(
    coverage_percent(data.frame(value=c(0, 1, 2, 0))), 50
))

bom <- list(
    bomFormat="CycloneDX", specVersion="1.5", version=1L,
    metadata=list(component=list(`bom-ref`="project")),
    components=list(list(`bom-ref`="dependency")),
    dependencies=list(
        list(ref="project", dependsOn=list("dependency")),
        list(ref="dependency", dependsOn=list())
    )
)
stopifnot(isTRUE(validate_cyclonedx_bom(bom)))
bad.bom <- bom
bad.bom$dependencies[[1L]]$dependsOn <- list("missing")
unknown.ref <- tryCatch({
    validate_cyclonedx_bom(bad.bom)
    FALSE
}, error=function(error) grepl("unknown", conditionMessage(error)))
stopifnot(isTRUE(unknown.ref))

page.calls <- 0L
pages <- collect_osv_pages(list(list(package="ape")), function(queries){
    page.calls <<- page.calls + 1L
    if(page.calls == 1L){
        list(list(vulns=list(list(id="A")), next_page_token="next"))
    } else{
        stopifnot(identical(queries[[1L]]$page_token, "next"))
        list(list(vulns=list(list(id="B"))))
    }
})
stopifnot(identical(
    vapply(pages[[1L]]$vulns, `[[`, character(1), "id"), c("A", "B")
))

# A generated file must match HEAD, even if it has already been staged.
check_generated_fixture <- function() {
    root <- tempfile("kfl1ou-generated-")
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
    git <- function(...) {
        status <- system2("git", c("-C", shQuote(root), vapply(
            c(...), shQuote, character(1)
        )), stdout = TRUE, stderr = TRUE)
        stopifnot(is.null(attr(status, "status")))
        invisible(status)
    }
    git("init", "--quiet")
    git("config", "user.name", "Generated-file fixture")
    git("config", "user.email", "fixture@example.invalid")
    git("config", "commit.gpgsign", "false")
    dir.create(file.path(root, "empty-hooks"))
    git("config", "core.hooksPath", file.path(root, "empty-hooks"))
    writeLines("export(original)", file.path(root, "NAMESPACE"))
    git("add", "NAMESPACE")
    git("commit", "--quiet", "-m", "Fixture baseline")
    stopifnot(isTRUE(check_generated(root)))
    rejects <- function() tryCatch({
        check_generated(root)
        FALSE
    }, error = function(error) grepl("differ from HEAD", conditionMessage(error)))
    writeLines("export(changed)", file.path(root, "NAMESPACE"))
    stopifnot(rejects())
    git("add", "NAMESPACE")
    stopifnot(rejects())
    git("commit", "--quiet", "-m", "Reviewed generated update")
    stopifnot(isTRUE(check_generated(root)))
    dir.create(file.path(root, "man"))
    writeLines("untracked generated page", file.path(root, "man", "new.Rd"))
    stopifnot(rejects())
    unlink(file.path(root, "man", "new.Rd"))
    writeLines("unrelated changes are allowed", file.path(root, "README.md"))
    stopifnot(isTRUE(check_generated(root)))
}
check_generated_fixture()

message("Maintenance-script fixtures passed.")
