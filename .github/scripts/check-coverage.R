coverage <- covr::package_coverage(quiet = FALSE)
coverage_data <- as.data.frame(coverage)

saveRDS(coverage, "coverage.rds")
utils::write.csv(coverage_data, "coverage.csv", row.names = FALSE)
covr::to_cobertura(coverage, filename = "coverage.xml")

coverage_percent <- function(data) {
  if (!nrow(data)) {
    return(NA_real_)
  }
  100 * mean(data$value > 0)
}

overall <- as.numeric(covr::percent_coverage(coverage))
file_floors <- c(
  "R/phylogeny_bootstrap.R" = 85,
  "R/univariate_sparse_backends.R" = 80,
  "R/rate_heterogeneity.R" = 80,
  "R/model_uncertainty.R" = 79,
  "R/measurement_error.R" = 74,
  "R/shift_configuration.R" = 80
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

format_percent <- function(value) {
  if (is.na(value)) "not applicable" else sprintf("%.2f%%", value)
}

summary_lines <- c(
  "# Coverage summary",
  "",
  sprintf("- Overall: **%s** (floor: %.1f%%)", format_percent(overall), overall_floor),
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
