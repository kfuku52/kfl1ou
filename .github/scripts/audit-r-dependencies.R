dependency.fields <- c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances")

description.field <- function(description, field) {
  if(field %in% names(description)) description[[field]] else NULL
}

parse.dependencies <- function(value) {
  if(is.null(value) || !length(value) || is.na(value) || !nzchar(value)){
    return(character())
  }
  entries <- trimws(unlist(strsplit(value, ",", fixed=TRUE)))
  packages <- trimws(sub("\\s*\\(.*$", "", entries))
  sort(unique(setdiff(packages[nzchar(packages)], "R")))
}

description.dependencies <- function(description, fields=dependency.fields) {
  sort(unique(unlist(lapply(fields, function(field) {
    parse.dependencies(description.field(description, field))
  }))))
}

project.description <- read.dcf("DESCRIPTION")[1, , drop=TRUE]
direct.dependencies <- description.dependencies(as.list(project.description))
pending <- direct.dependencies
records <- list()
dependency.edges <- list()
base.packages <- character()
missing.packages <- character()

while(length(pending)){
  package <- pending[[1L]]
  pending <- pending[-1L]
  if(package %in% c(names(records), base.packages, missing.packages)){
    next
  }

  description <- tryCatch(
    suppressWarnings(utils::packageDescription(package)),
    error=function(error) NULL
  )
  if(is.null(description) ||
     is.null(description.field(description, "Package")) ||
     is.na(description.field(description, "Package"))){
    missing.packages <- c(missing.packages, package)
    next
  }

  dependencies <- description.dependencies(
    description,
    c("Depends", "Imports", "LinkingTo")
  )
  if(identical(description.field(description, "Priority"), "base")){
    base.packages <- c(base.packages, package)
    next
  }

  records[[package]] <- list(
    name=package,
    version=as.character(description.field(description, "Version")),
    license=if(is.null(description.field(description, "License")))
      NA_character_ else as.character(description.field(description, "License")),
    direct=package %in% direct.dependencies
  )
  dependency.edges[[package]] <- dependencies
  pending <- sort(unique(c(
    pending,
    setdiff(
      dependencies,
      c(names(records), base.packages, missing.packages)
    )
  )))
}

if(length(missing.packages)){
  stop(
    "Dependency audit could not inspect installed packages: ",
    paste(sort(unique(missing.packages)), collapse=", "),
    call.=FALSE
  )
}

records <- records[sort(names(records))]
project.ref <- paste0(
  "pkg:cran/",
  unname(project.description[["Package"]]),
  "@",
  unname(project.description[["Version"]])
)
component.refs <- stats::setNames(
  vapply(records, function(record) {
    paste0(
      "pkg:cran/",
      utils::URLencode(record$name, reserved=TRUE),
      "@",
      utils::URLencode(record$version, reserved=TRUE)
    )
  }, character(1)),
  names(records)
)

components <- lapply(names(records), function(package) {
  record <- records[[package]]
  list(
    type="library",
    `bom-ref`=unname(component.refs[[package]]),
    name=record$name,
    version=record$version,
    purl=unname(component.refs[[package]]),
    licenses=if(is.na(record$license)) NULL else list(list(
      license=list(name=record$license)
    )),
    properties=list(list(
      name="kfl1ou:direct-dependency",
      value=if(record$direct) "true" else "false"
    ))
  )
})

dependencies <- c(list(list(
  ref=project.ref,
  dependsOn=I(unname(component.refs[
    intersect(direct.dependencies, names(records))
  ]))
)), lapply(names(records), function(package) {
  children <- intersect(dependency.edges[[package]], names(records))
  list(
    ref=unname(component.refs[[package]]),
    dependsOn=I(unname(component.refs[children]))
  )
}))

bom <- list(
  bomFormat="CycloneDX",
  specVersion="1.5",
  version=1,
  metadata=list(
    timestamp=format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC"),
    component=list(
      type="application",
      name=unname(project.description[["Package"]]),
      version=unname(project.description[["Version"]]),
      `bom-ref`=project.ref
    )
  ),
  components=components,
  dependencies=dependencies
)
jsonlite::write_json(
  bom,
  "r-dependency-sbom.cdx.json",
  auto_unbox=TRUE,
  pretty=TRUE,
  null="null"
)
serialized.bom <- jsonlite::fromJSON(
  "r-dependency-sbom.cdx.json",
  simplifyVector=FALSE
)
invalid.dependency.arrays <- vapply(
  serialized.bom$dependencies,
  function(dependency) !is.list(dependency$dependsOn),
  logical(1)
)
if(any(invalid.dependency.arrays)){
  stop("CycloneDX dependsOn values must serialize as JSON arrays.", call.=FALSE)
}

queries <- lapply(records, function(record) {
  list(
    package=list(ecosystem="CRAN", name=record$name),
    version=record$version
  )
})
query.osv.batch <- function(batch.queries) {
  payload <- jsonlite::toJSON(
    list(queries=unname(batch.queries)),
    auto_unbox=TRUE,
    null="null"
  )
  handle <- curl::new_handle(
    postfields=charToRaw(payload),
    httpheader=c("Content-Type"="application/json"),
    timeout=60
  )
  response <- curl::curl_fetch_memory(
    "https://api.osv.dev/v1/querybatch",
    handle=handle
  )
  if(response$status_code != 200L){
    stop(
      "OSV query failed with HTTP status ", response$status_code,
      call.=FALSE
    )
  }
  result <- jsonlite::fromJSON(
    rawToChar(response$content),
    simplifyVector=FALSE
  )
  if(length(result$results) != length(batch.queries)){
    stop("OSV returned an incomplete batch response.", call.=FALSE)
  }
  result$results
}

osv.results <- rep(list(list(vulns=list())), length(queries))
pending.indices <- seq_along(queries)
pending.queries <- queries
seen.tokens <- rep(list(character()), length(queries))
page.round <- 0L
while(length(pending.indices)){
  page.round <- page.round + 1L
  if(page.round > 100L){
    stop("OSV pagination exceeded 100 rounds.", call.=FALSE)
  }
  page.results <- query.osv.batch(pending.queries)
  next.indices <- integer()
  next.queries <- list()
  for(batch.index in seq_along(page.results)){
    record.index <- pending.indices[[batch.index]]
    page <- page.results[[batch.index]]
    page.vulnerabilities <- if(is.null(page$vulns)) list() else page$vulns
    osv.results[[record.index]]$vulns <- c(
      osv.results[[record.index]]$vulns,
      page.vulnerabilities
    )
    token <- page$next_page_token
    if(!is.null(token) && length(token) == 1L && nzchar(token)){
      if(token %in% seen.tokens[[record.index]]){
        stop("OSV returned a repeated pagination token.", call.=FALSE)
      }
      seen.tokens[[record.index]] <- c(
        seen.tokens[[record.index]], token
      )
      next.indices <- c(next.indices, record.index)
      next.query <- queries[[record.index]]
      next.query$page_token <- token
      next.queries[[length(next.queries) + 1L]] <- next.query
    }
  }
  pending.indices <- next.indices
  pending.queries <- next.queries
}

vulnerability.ids <- sort(unique(unlist(lapply(osv.results, function(result) {
  vulnerabilities <- result$vulns
  if(is.null(vulnerabilities) || !length(vulnerabilities)){
    return(character())
  }
  vapply(vulnerabilities, function(vulnerability) {
    if(is.null(vulnerability$id) || length(vulnerability$id) != 1L){
      stop("OSV batch response contained a vulnerability without one ID.",
           call.=FALSE)
    }
    as.character(vulnerability$id)
  }, character(1))
}), use.names=FALSE)))

fetch.osv.vulnerability <- function(id, attempts=3L) {
  url <- paste0(
    "https://api.osv.dev/v1/vulns/",
    curl::curl_escape(id)
  )
  last.error <- NULL
  for(attempt in seq_len(attempts)){
    response <- tryCatch(
      curl::curl_fetch_memory(
        url,
        handle=curl::new_handle(timeout=60)
      ),
      error=function(error) {
        last.error <<- conditionMessage(error)
        NULL
      }
    )
    if(!is.null(response) && response$status_code == 200L){
      record <- jsonlite::fromJSON(
        rawToChar(response$content),
        simplifyVector=FALSE
      )
      if(is.null(record$id) || !identical(as.character(record$id), id)){
        stop("OSV vulnerability detail response had a mismatched ID.",
             call.=FALSE)
      }
      return(record)
    }
    if(!is.null(response)){
      last.error <- paste("HTTP status", response$status_code)
    }
    if(attempt < attempts){
      Sys.sleep(attempt)
    }
  }
  stop(
    "OSV vulnerability detail request failed for ", id, ": ",
    if(is.null(last.error)) "unknown error" else last.error,
    call.=FALSE
  )
}

vulnerability.details <- stats::setNames(
  lapply(vulnerability.ids, fetch.osv.vulnerability),
  vulnerability.ids
)

vulnerability.rows <- list()
row.index <- 0L
for(index in seq_along(records)){
  vulnerabilities <- osv.results[[index]]$vulns
  if(is.null(vulnerabilities) || !length(vulnerabilities)){
    next
  }
  ids <- vapply(vulnerabilities, function(vulnerability) {
    as.character(vulnerability$id)
  }, character(1))
  batch.modified <- vapply(vulnerabilities, function(vulnerability) {
    if(is.null(vulnerability$modified)) "" else
      as.character(vulnerability$modified)
  }, character(1))
  keep <- !duplicated(ids)
  ids <- ids[keep]
  batch.modified <- batch.modified[keep]
  for(vulnerability.index in seq_along(ids)){
    vulnerability <- vulnerability.details[[ids[[vulnerability.index]]]]
    if(!is.null(vulnerability$withdrawn) &&
       length(vulnerability$withdrawn) == 1L &&
       nzchar(as.character(vulnerability$withdrawn))){
      next
    }
    row.index <- row.index + 1L
    severity <- if(is.null(vulnerability$severity)) character() else
      vapply(vulnerability$severity, function(item) {
        if(is.null(item$score)) "" else as.character(item$score)
      }, character(1))
    vulnerability.rows[[row.index]] <- data.frame(
      package=records[[index]]$name,
      version=records[[index]]$version,
      id=as.character(vulnerability$id),
      aliases=paste(unlist(vulnerability$aliases), collapse=";"),
      severity=paste(severity[nzchar(severity)], collapse=";"),
      summary=if(is.null(vulnerability$summary)) "" else
        as.character(vulnerability$summary),
      modified=if(is.null(vulnerability$modified))
        batch.modified[[vulnerability.index]] else
          as.character(vulnerability$modified),
      stringsAsFactors=FALSE
    )
  }
}
vulnerabilities <- if(length(vulnerability.rows)){
  do.call(rbind, vulnerability.rows)
} else{
  data.frame(
    package=character(),
    version=character(),
    id=character(),
    aliases=character(),
    severity=character(),
    summary=character(),
    modified=character(),
    stringsAsFactors=FALSE
  )
}
utils::write.csv(
  vulnerabilities,
  "r-dependency-vulnerabilities.csv",
  row.names=FALSE
)

summary.lines <- c(
  "# R dependency audit",
  "",
  paste0("- Direct dependencies: **", length(direct.dependencies), "**"),
  paste0("- Installed dependency closure in SBOM: **", length(records), "**"),
  paste0("- Active OSV findings: **", nrow(vulnerabilities), "**"),
  paste0(
    "- Generated: ",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz="UTC")
  )
)
writeLines(summary.lines, "r-dependency-audit.md")
step.summary <- Sys.getenv("GITHUB_STEP_SUMMARY")
if(nzchar(step.summary)){
  cat(
    paste(summary.lines, collapse="\n"),
    "\n",
    file=step.summary,
    append=TRUE
  )
}
cat(paste(summary.lines, collapse="\n"), "\n")

if(nrow(vulnerabilities)){
  stop(
    "OSV reported active vulnerabilities in the installed R dependency closure.",
    call.=FALSE
  )
}
