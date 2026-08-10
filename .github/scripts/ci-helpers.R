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

coverage_percent <- function(data) {
  if(!nrow(data)) return(NA_real_)
  100 * mean(data$value > 0)
}

format_percent <- function(value) {
  if(is.na(value)) "not applicable" else sprintf("%.2f%%", value)
}

validate_cyclonedx_bom <- function(bom) {
  required <- c("bomFormat", "specVersion", "version", "metadata",
                "components", "dependencies")
  if(!all(required %in% names(bom)) || !identical(bom$bomFormat, "CycloneDX")){
    stop("invalid CycloneDX document structure.", call.=FALSE)
  }
  project.ref <- bom$metadata$component[["bom-ref"]]
  component.refs <- vapply(bom$components, `[[`, character(1), "bom-ref")
  if(anyDuplicated(component.refs) || is.null(project.ref) || !nzchar(project.ref)){
    stop("CycloneDX bom-ref values must be non-empty and unique.", call.=FALSE)
  }
  known.refs <- c(project.ref, component.refs)
  dependency.refs <- vapply(bom$dependencies, `[[`, character(1), "ref")
  if(anyDuplicated(dependency.refs) || !setequal(dependency.refs, known.refs)){
    stop("CycloneDX dependency refs must cover every component exactly once.",
         call.=FALSE)
  }
  for(dependency in bom$dependencies){
    if(!is.list(dependency$dependsOn)){
      stop("CycloneDX dependsOn values must serialize as JSON arrays.", call.=FALSE)
    }
    children <- unlist(dependency$dependsOn, use.names=FALSE)
    if(any(!children %in% known.refs)){
      stop("CycloneDX dependency graph contains an unknown bom-ref.", call.=FALSE)
    }
  }
  invisible(TRUE)
}

collect_osv_pages <- function(queries, fetch, max.rounds=100L) {
  results <- rep(list(list(vulns=list())), length(queries))
  pending.indices <- seq_along(queries)
  pending.queries <- queries
  seen.tokens <- rep(list(character()), length(queries))
  page.round <- 0L
  while(length(pending.indices)){
    page.round <- page.round + 1L
    if(page.round > max.rounds){
      stop("OSV pagination exceeded the configured round limit.", call.=FALSE)
    }
    pages <- fetch(pending.queries)
    if(length(pages) != length(pending.queries)){
      stop("OSV returned an incomplete batch response.", call.=FALSE)
    }
    next.indices <- integer()
    next.queries <- list()
    for(batch.index in seq_along(pages)){
      record.index <- pending.indices[[batch.index]]
      page <- pages[[batch.index]]
      vulnerabilities <- if(is.null(page$vulns)) list() else page$vulns
      results[[record.index]]$vulns <- c(
        results[[record.index]]$vulns, vulnerabilities
      )
      token <- page$next_page_token
      if(!is.null(token) && length(token) == 1L && nzchar(token)){
        if(token %in% seen.tokens[[record.index]]){
          stop("OSV returned a repeated pagination token.", call.=FALSE)
        }
        seen.tokens[[record.index]] <- c(seen.tokens[[record.index]], token)
        next.indices <- c(next.indices, record.index)
        next.query <- queries[[record.index]]
        next.query$page_token <- token
        next.queries[[length(next.queries) + 1L]] <- next.query
      }
    }
    pending.indices <- next.indices
    pending.queries <- next.queries
  }
  results
}
