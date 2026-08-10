library(ape)
library(kfl1ou)

benchmark_profile <- Sys.getenv("KFL1OU_BENCHMARK_PROFILE", "smoke")
if (!benchmark_profile %in% c("smoke", "extended")) {
  stop("KFL1OU_BENCHMARK_PROFILE must be `smoke` or `extended`.", call. = FALSE)
}
max_regression_factor <- as.numeric(
  Sys.getenv("KFL1OU_BENCHMARK_MAX_FACTOR", "3")
)
if (!is.finite(max_regression_factor) || max_regression_factor <= 1) {
  stop("KFL1OU_BENCHMARK_MAX_FACTOR must be a finite number greater than one.",
       call. = FALSE)
}

benchmark_budget <- function(name, median_seconds, baseline_seconds,
                             absolute_seconds) {
  relative_limit <- baseline_seconds * max_regression_factor
  effective_limit <- min(relative_limit, absolute_seconds)
  if (median_seconds > effective_limit) {
    stop(
      sprintf(
        paste0(
          "%s median runtime %.3fs exceeds the %.3fs effective budget ",
          "(baseline %.3fs x %.2f, absolute ceiling %.3fs)."
        ),
        name, median_seconds, effective_limit, baseline_seconds,
        max_regression_factor, absolute_seconds
      ),
      call. = FALSE
    )
  }
  list(
    baseline_seconds = baseline_seconds,
    regression_ratio = median_seconds / baseline_seconds,
    effective_limit_seconds = effective_limit,
    absolute_limit_seconds = absolute_seconds
  )
}

make_data <- function(n_tips, seed) {
  set.seed(seed)
  tree <- reorder.phylo(rcoal(n_tips), "postorder")
  trait <- rTraitCont(tree, model = "OU", alpha = 1, sigma = 0.4)
  names(trait) <- tree$tip.label
  list(tree = tree, trait = trait)
}

benchmark_case <- function(name, data, strategy, max_shifts,
                           baseline_seconds, absolute_seconds,
                           iterations = 3L) {
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
  budget <- benchmark_budget(
    name, median_seconds, baseline_seconds, absolute_seconds
  )

  data.frame(
    benchmark = name,
    strategy = strategy,
    tips = length(data$tree$tip.label),
    max_shifts = max_shifts,
    iterations = iterations,
    median_seconds = median_seconds,
    baseline_seconds = budget$baseline_seconds,
    regression_ratio = budget$regression_ratio,
    effective_limit_seconds = budget$effective_limit_seconds,
    absolute_limit_seconds = budget$absolute_limit_seconds,
    min_seconds = min(elapsed),
    max_seconds = max(elapsed),
    score = reference$score,
    selected_shifts = reference$nShifts,
    stringsAsFactors = FALSE
  )
}

benchmark_full_covariance_case <- function(baseline_seconds,
                                           absolute_seconds,
                                           iterations = 3L) {
  set.seed(3L)
  tree <- reorder.phylo(rcoal(84L), "postorder")
  Y <- matrix(
    rnorm(84L * 3L), nrow = 84L, ncol = 3L,
    dimnames = list(tree$tip.label, paste0("trait", seq_len(3L)))
  )
  run_once <- function() {
    suppressWarnings(
      fit_OU(
        tree,
        Y,
        integer(),
        criterion = "BIC",
        trait.covariance = "full",
        alpha.structure = "diagonal",
        alpha.lower = 0.7,
        alpha.upper = 0.7,
        likelihood.engine = "auto",
        optimizer.starts = 1L,
        compute.hessian = FALSE
      )
    )
  }

  reference <- run_once()
  stopifnot(
    identical(reference$likelihood.engine, "dense"),
    identical(
      reference$diagnostics$engine.selection$reason,
      "estimated-dense-cost"
    )
  )
  elapsed <- numeric(iterations)
  for (index in seq_len(iterations)) {
    timing <- system.time(candidate <- run_once())
    elapsed[[index]] <- unname(timing[["elapsed"]])
    stopifnot(
      identical(candidate$likelihood.engine, reference$likelihood.engine),
      isTRUE(all.equal(
        candidate$joint.logLik,
        reference$joint.logLik,
        tolerance = 1e-10
      )),
      isTRUE(all.equal(candidate$score, reference$score, tolerance = 1e-10))
    )
  }

  median_seconds <- stats::median(elapsed)
  budget <- benchmark_budget(
    "automatic-full-covariance", median_seconds, baseline_seconds,
    absolute_seconds
  )

  data.frame(
    benchmark = "automatic-full-covariance",
    strategy = reference$likelihood.engine,
    tips = length(tree$tip.label),
    max_shifts = 0L,
    iterations = iterations,
    median_seconds = median_seconds,
    baseline_seconds = budget$baseline_seconds,
    regression_ratio = budget$regression_ratio,
    effective_limit_seconds = budget$effective_limit_seconds,
    absolute_limit_seconds = budget$absolute_limit_seconds,
    min_seconds = min(elapsed),
    max_seconds = max(elapsed),
    score = reference$score,
    selected_shifts = reference$nShifts,
    stringsAsFactors = FALSE
  )
}

benchmark_pruning_cache_case <- function(baseline_seconds,
                                         absolute_seconds,
                                         iterations = 2L) {
  set.seed(4L)
  tree <- reorder.phylo(rcoal(300L), "postorder")
  Y <- matrix(
    rnorm(600L), nrow=300L, ncol=2L,
    dimnames=list(tree$tip.label, c("trait1", "trait2"))
  )
  run_once <- function() {
    suppressWarnings(fit_OU(
      tree, Y, integer(), criterion="BIC",
      trait.covariance="full", alpha.structure="diagonal",
      alpha.lower=c(0.6, 0.9), alpha.upper=c(0.6, 0.9),
      likelihood.engine="pruning", optimizer.starts=1L,
      compute.hessian=FALSE
    ))
  }
  reference <- run_once()
  stopifnot(
    identical(reference$likelihood.engine, "pruning"),
    identical(reference$diagnostics$tree.cache$dense, FALSE)
  )
  elapsed <- numeric(iterations)
  for(index in seq_len(iterations)) {
    timing <- system.time(candidate <- run_once())
    elapsed[[index]] <- unname(timing[["elapsed"]])
    stopifnot(isTRUE(all.equal(
      candidate$joint.logLik, reference$joint.logLik, tolerance=1e-10
    )))
  }
  median_seconds <- stats::median(elapsed)
  budget <- benchmark_budget(
    "large-explicit-pruning", median_seconds, baseline_seconds,
    absolute_seconds
  )
  data.frame(
    benchmark="large-explicit-pruning",
    strategy="pruning",
    tips=length(tree$tip.label),
    max_shifts=0L,
    iterations=iterations,
    median_seconds=median_seconds,
    baseline_seconds=budget$baseline_seconds,
    regression_ratio=budget$regression_ratio,
    effective_limit_seconds=budget$effective_limit_seconds,
    absolute_limit_seconds=budget$absolute_limit_seconds,
    min_seconds=min(elapsed),
    max_seconds=max(elapsed),
    score=reference$score,
    selected_shifts=reference$nShifts,
    stringsAsFactors=FALSE
  )
}

results <- list(
  benchmark_case(
    "certified-small-search",
    make_data(12L, 1L),
    "exhaustive",
    2L,
    0.5,
    2
  ),
  benchmark_case(
    "ensemble-medium-search",
    make_data(30L, 2L),
    "ensemble",
    3L,
    0.15,
    1
  ),
  benchmark_full_covariance_case(1.5, 6)
)
if(identical(benchmark_profile, "extended")) {
  results[[length(results) + 1L]] <- benchmark_pruning_cache_case(2.5, 10)
}
results <- do.call(rbind, results)

utils::write.csv(results, "benchmark-results.csv", row.names = FALSE)
print(results, row.names = FALSE)

step_summary <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(step_summary)) {
  lines <- c(
    "# Search benchmark summary",
    "",
    sprintf("Profile: **%s**; maximum baseline regression: **%.2fx**.",
            benchmark_profile, max_regression_factor),
    "",
    "| Benchmark | Strategy | Tips | Median | Baseline | Ratio | Limit | Selected shifts |",
    "|---|---|---:|---:|---:|---:|---:|---:|",
    sprintf(
      "| %s | %s | %d | %.3f | %.3f | %.2fx | %.3f | %d |",
      results$benchmark,
      results$strategy,
      results$tips,
      results$median_seconds,
      results$baseline_seconds,
      results$regression_ratio,
      results$effective_limit_seconds,
      results$selected_shifts
    ),
    ""
  )
  cat(paste(lines, collapse = "\n"), "\n", file = step_summary, append = TRUE)
}
