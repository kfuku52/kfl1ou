library(ape)
library(kfl1ou)

make_data <- function(n_tips, seed) {
  set.seed(seed)
  tree <- reorder.phylo(rcoal(n_tips), "postorder")
  trait <- rTraitCont(tree, model = "OU", alpha = 1, sigma = 0.4)
  names(trait) <- tree$tip.label
  list(tree = tree, trait = trait)
}

benchmark_case <- function(name, data, strategy, max_shifts,
                           max_median_seconds, iterations = 3L) {
  run_once <- function() {
    set.seed(20260723)
    suppressWarnings(
      estimate_shift_configuration(
        data$tree,
        data$trait,
        max.nShifts = max_shifts,
        search.strategy = strategy,
        ensemble.replicates = if (strategy == "ensemble") 3L else 8L,
        quietly = TRUE
      )
    )
  }

  reference <- run_once()
  elapsed <- numeric(iterations)
  for (index in seq_len(iterations)) {
    timing <- system.time(candidate <- run_once())
    elapsed[[index]] <- unname(timing[["elapsed"]])
    stopifnot(
      identical(candidate$shift.configuration, reference$shift.configuration),
      isTRUE(all.equal(candidate$score, reference$score, tolerance = 1e-10)),
      identical(candidate$search.diagnostics$strategy, strategy)
    )
  }

  median_seconds <- stats::median(elapsed)
  if (median_seconds > max_median_seconds) {
    stop(
      sprintf(
        "%s median runtime %.3fs exceeds the %.3fs regression budget.",
        name,
        median_seconds,
        max_median_seconds
      ),
      call. = FALSE
    )
  }

  data.frame(
    benchmark = name,
    strategy = strategy,
    tips = length(data$tree$tip.label),
    max_shifts = max_shifts,
    iterations = iterations,
    median_seconds = median_seconds,
    max_median_seconds = max_median_seconds,
    min_seconds = min(elapsed),
    max_seconds = max(elapsed),
    score = reference$score,
    selected_shifts = reference$nShifts,
    stringsAsFactors = FALSE
  )
}

results <- rbind(
  benchmark_case(
    "certified-small-search",
    make_data(12L, 1L),
    "exhaustive",
    2L,
    5
  ),
  benchmark_case(
    "ensemble-medium-search",
    make_data(30L, 2L),
    "ensemble",
    3L,
    2
  )
)

utils::write.csv(results, "benchmark-results.csv", row.names = FALSE)
print(results, row.names = FALSE)

step_summary <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(step_summary)) {
  lines <- c(
    "# Search benchmark summary",
    "",
    "| Benchmark | Strategy | Tips | Median seconds | Budget seconds | Selected shifts |",
    "|---|---|---:|---:|---:|---:|",
    sprintf(
      "| %s | %s | %d | %.3f | %.3f | %d |",
      results$benchmark,
      results$strategy,
      results$tips,
      results$median_seconds,
      results$max_median_seconds,
      results$selected_shifts
    ),
    ""
  )
  cat(paste(lines, collapse = "\n"), "\n", file = step_summary, append = TRUE)
}
