source(file.path(".github", "scripts", "check-metadata.R"), local = TRUE)
source(file.path(".github", "scripts", "ci-helpers.R"), local = TRUE)

fixture <- tempfile("kfl1ou-metadata-")
dir.create(fixture)
on.exit(unlink(fixture, recursive = TRUE, force = TRUE), add = TRUE)

writeLines(c(
    "Package: kfl1ou",
    "Version: 9.8.7",
    paste("Date:", Sys.Date())
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

message("Maintenance-script fixtures passed.")
