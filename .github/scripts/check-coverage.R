source(file.path(".github", "scripts", "ci-helpers.R"), local=TRUE)

coverage <- covr::package_coverage(quiet = FALSE)
coverage_data <- as.data.frame(coverage)

saveRDS(coverage, "coverage.rds")
utils::write.csv(coverage_data, "coverage.csv", row.names = FALSE)
covr::to_cobertura(coverage, filename = "coverage.xml")

overall <- as.numeric(covr::percent_coverage(coverage))
namespace <- trimws(readLines("NAMESPACE", warn=FALSE))
export.lines <- grep("^export\\([^)]+\\)$", namespace, value=TRUE)
s3.lines <- grep("^S3method\\([^,]+,[^)]+\\)$", namespace, value=TRUE)
public.api <- sort(unique(c(
  sub("^export\\(([^)]+)\\)$", "\\1", export.lines),
  sub("^S3method\\(([^,]+),([^)]+)\\)$", "\\1.\\2", s3.lines)
)))
coverage.functions <- as.character(coverage_data$functions)
public.api.results <- data.frame(
  function_name=public.api,
  expressions=vapply(
    public.api,
    function(name) sum(coverage.functions == name, na.rm=TRUE),
    integer(1)
  ),
  covered_expressions=vapply(
    public.api,
    function(name) sum(
      coverage.functions == name & coverage_data$value > 0,
      na.rm=TRUE
    ),
    integer(1)
  ),
  stringsAsFactors=FALSE
)
public.api.results$coverage <- ifelse(
  public.api.results$expressions == 0L,
  NA_real_,
  100 * public.api.results$covered_expressions /
    public.api.results$expressions
)
utils::write.csv(
  public.api.results,
  "public-api-coverage.csv",
  row.names=FALSE
)
zero.public.api <- public.api.results$function_name[
  public.api.results$expressions == 0L |
    public.api.results$covered_expressions == 0L
]
file_floors <- c(
  "R/phylogeny_bootstrap.R" = 85,
  "R/univariate_sparse_backends.R" = 80,
  "R/rate_heterogeneity.R" = 80,
  "R/model_uncertainty.R" = 79,
  "R/measurement_error.R" = 74,
  "R/shift_configuration.R" = 80,
  "R/multivariate_covariance.R" = 85,
  "R/pruning_likelihood.R" = 75,
  "R/convergent_regions.R" = 85,
  "R/model_api.R" = 85
)
file_results <- data.frame(
  file = names(file_floors),
  coverage = vapply(
    names(file_floors),
    function(path) coverage_percent(
      coverage_data[coverage_data$filename == path, , drop = FALSE]
    ),
    numeric(1)
  ),
  floor = unname(file_floors),
  stringsAsFactors = FALSE
)

changed_coverage <- function(data, base_sha) {
  empty <- list(percent = NA_real_, expressions = 0L, files = character())
  if (!nzchar(base_sha) || grepl("^0+$", base_sha)) {
    return(empty)
  }

  exists <- suppressWarnings(
    system2(
      "git", c("cat-file", "-e", paste0(base_sha, "^{commit}")),
      stdout = FALSE, stderr = FALSE
    )
  )
  if (!identical(exists, 0L)) {
    return(empty)
  }

  diff <- system2(
    "git",
    c("diff", "--unified=0", "--no-color", paste0(base_sha, "..HEAD"), "--", "R"),
    stdout = TRUE,
    stderr = TRUE
  )
  current_file <- NA_character_
  changed <- list()

  for (line in diff) {
    if (startsWith(line, "+++ b/")) {
      current_file <- sub("^\\+\\+\\+ b/", "", line)
      next
    }
    if (!startsWith(line, "@@") || is.na(current_file)) {
      next
    }

    match <- regexec(
      "@@ -[0-9]+(?:,[0-9]+)? \\+([0-9]+)(?:,([0-9]+))? @@",
      line,
      perl = TRUE
    )
    parts <- regmatches(line, match)[[1L]]
    if (!length(parts)) {
      next
    }
    first <- as.integer(parts[[2L]])
    count <- if (length(parts) >= 3L && nzchar(parts[[3L]])) {
      as.integer(parts[[3L]])
    } else {
      1L
    }
    if (count > 0L) {
      changed[[current_file]] <- unique(c(
        changed[[current_file]],
        seq.int(first, length.out = count)
      ))
    }
  }

  covered_rows <- integer()
  for (file in intersect(names(changed), unique(data$filename))) {
    candidates <- which(data$filename == file)
    hit <- vapply(
      candidates,
      function(index) {
        any(changed[[file]] >= data$first_line[[index]] &
          changed[[file]] <= data$last_line[[index]])
      },
      logical(1)
    )
    covered_rows <- c(covered_rows, candidates[hit])
  }
  covered_rows <- unique(covered_rows)
  if (!length(covered_rows)) {
    return(empty)
  }

  list(
    percent = coverage_percent(data[covered_rows, , drop = FALSE]),
    expressions = length(covered_rows),
    files = intersect(names(changed), unique(data$filename))
  )
}

patch <- changed_coverage(coverage_data, Sys.getenv("BASE_SHA"))
overall_floor <- 80.5
patch_floor <- 80

summary_lines <- c(
  "# Coverage summary",
  "",
  sprintf("- Overall: **%s** (floor: %.1f%%)", format_percent(overall), overall_floor),
  if(length(zero.public.api)){
    paste0(
      "- Public API functions with zero covered expressions: **",
      paste(zero.public.api, collapse=", "),
      "**"
    )
  } else{
    sprintf(
      "- Public API functions with zero covered expressions: **none** (%d checked)",
      nrow(public.api.results)
    )
  },
  if (is.na(patch$percent)) {
    "- Changed executable expressions: **not applicable**"
  } else {
    sprintf(
      "- Changed executable expressions: **%s** across %d expressions (floor: %.0f%%)",
      format_percent(patch$percent), patch$expressions, patch_floor
    )
  },
  "",
  "| Critical file | Coverage | Floor |",
  "|---|---:|---:|",
  sprintf(
    "| `%s` | %s | %.0f%% |",
    file_results$file,
    vapply(file_results$coverage, format_percent, character(1)),
    file_results$floor
  ),
  ""
)
writeLines(summary_lines, "coverage-summary.md")

step_summary <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(step_summary)) {
  cat(paste(summary_lines, collapse = "\n"), "\n", file = step_summary, append = TRUE)
}
cat(paste(summary_lines, collapse = "\n"), "\n")

failures <- character()
if(length(zero.public.api)){
  failures <- c(
    failures,
    paste0(
      "Public API functions must have at least one covered expression: ",
      paste(zero.public.api, collapse=", ")
    )
  )
}
if (overall < overall_floor) {
  failures <- c(
    failures,
    sprintf("Overall coverage %.2f%% is below %.1f%%", overall, overall_floor)
  )
}
below_file_floor <- is.na(file_results$coverage) |
  file_results$coverage < file_results$floor
if (any(below_file_floor)) {
  failures <- c(
    failures,
    sprintf(
      "%s coverage %s is below %.0f%%",
      file_results$file[below_file_floor],
      vapply(
        file_results$coverage[below_file_floor],
        format_percent,
        character(1)
      ),
      file_results$floor[below_file_floor]
    )
  )
}
if (!is.na(patch$percent) && patch$percent < patch_floor) {
  failures <- c(
    failures,
    sprintf(
      "Changed-expression coverage %.2f%% is below %.0f%%",
      patch$percent, patch_floor
    )
  )
}

if (length(failures)) {
  stop(paste(failures, collapse = "\n"), call. = FALSE)
}
