test_that("tree ensembles align unnamed errors with the original response rows", {
  dat <- small_lizard_data(n_tips=8, traits=1:2)
  Y <- dat$Y[8:1, , drop=FALSE]
  errors <- matrix(seq(.01, .8, length.out=16), 8, 2)
  named <- errors
  dimnames(named) <- dimnames(Y)
  trees <- list(dat$tree, dat$tree)
  trees[[2]]$tip.label <- rev(trees[[2]]$tip.label)
  for(input in list(errors, errors[, 1],
                   setNames(as.data.frame(errors), colnames(Y)))){
    expected.error <- if(is.null(dim(input))) named[, 1] else named
    if(is.null(dim(input))) names(expected.error) <- rownames(Y)
    expected <- fit_l1ou_tree_ensemble(
      trees, Y, input_error=expected.error, max.nShifts=0, criterion="BIC",
      alpha.lower=.7, alpha.upper=.7
    )
    actual <- fit_l1ou_tree_ensemble(
      trees, Y, input_error=input, max.nShifts=0, criterion="BIC",
      alpha.lower=.7, alpha.upper=.7
    )
    for(i in 1:2){
      expect_equal(actual$fits[[i]]$l1ou.options$input_error,
                   expected$fits[[i]]$l1ou.options$input_error)
      expect_equal(as.numeric(logLik(actual$fits[[i]])),
                   as.numeric(logLik(expected$fits[[i]])))
    }
  }
  options <- expected$fits[[1]]$l1ou.options
  options$input_error <- errors
  nested <- fit_l1ou_tree_ensemble(trees, Y, l1ou.options=options, nCores=2L)
  reference <- fit_l1ou_tree_ensemble(
    trees, Y, input_error=named, max.nShifts=0, criterion="BIC",
    alpha.lower=.7, alpha.upper=.7
  )
  expect_equal(lapply(nested$fits, function(x) x$l1ou.options$input_error),
               lapply(reference$fits, function(x) x$l1ou.options$input_error))
})

test_that("model averaging excludes non-finite scores before choosing candidates", {
  dat <- small_lizard_data(n_tips=8)
  fit <- fit_OU(dat$tree, dat$Y, integer(), criterion="BIC",
                alpha.lower=.7, alpha.upper=.7)
  fit$profile <- list(scores=c(-Inf, NA, Inf, fit$score),
                      configurations=list(integer(), 1L, 2L, integer()))
  average <- model_average_l1ou(fit)
  expect_equal(average$weights, 1)
  expect_equal(average$configurations, list(integer()))
  expect_equal(average$mu, fit$mu)
  fit$profile <- list(scores=c(1, 2), configurations=list(integer()))
  expect_error(model_average_l1ou(fit), "one score per configuration")
  fit$profile <- list(scores=1, configurations=list(1.5))
  expect_error(model_average_l1ou(fit), "integer edge indices")
})

test_that("configuration lookup allows gaps in evaluated shift counts", {
  model <- structure(list(profile=list(
    scores=c(10, 2, 1),
    configurations=list(integer(), 1:5, 2:6)
  )), class="l1ou")
  expect_equal(get_shift_configuration(model, 5L), 2:6)
  expect_equal(get_shift_configuration(model, 0L), integer())
  expect_error(get_shift_configuration(model, 3L), "There is no configuration")
})

test_that("uncertainty summaries reject invalid edges before integer conversion", {
  dat <- small_lizard_data(n_tips=8)
  for(tree in list(NULL, dat$tree)){
    for(shifts in list(1.5, "1", TRUE, NA_real_, Inf, -1, c(1L, 1L))){
      expect_error(summarize_shift_uncertainty(
        list(all.shifts=list(shifts), detection.rate=numeric(14)), tree
      ), "invalid shift-edge")
    }
  }
  valid <- summarize_shift_uncertainty(
    list(all.shifts=list(integer(), c(1L, 3L)), detection.rate=numeric(14))
  )
  expect_equal(unname(valid$edge.inclusion[c(1, 3)]), c(.5, .5))
  empty <- summarize_shift_uncertainty(list(all.shifts=list(integer())))
  expect_equal(empty$successful, 1L)
  expect_length(empty$edge.inclusion, 0L)
})

test_that("placeholder tables accept positional data without losing explicit labels", {
  dat <- small_lizard_data(n_tips=8, traits=1:2)
  Y <- unname(dat$Y)
  for(input in list(Y, as.data.frame(Y))){
    table <- kfl1ou:::as_placeholder_trait_table(dat$tree, input)
    expect_equal(rownames(table), dat$tree$tip.label)
    expect_equal(unname(as.matrix(table)), Y)
  }
  named <- dat$Y[8:1, , drop=FALSE]
  expect_equal(as.matrix(kfl1ou:::as_placeholder_trait_table(dat$tree, named)), dat$Y)
  rownames(named)[1] <- "unknown"
  expect_error(kfl1ou:::as_placeholder_trait_table(dat$tree, named), "rownames")
  directory <- tempfile("kfl1ou-positional-output-")
  on.exit(unlink(directory, recursive=TRUE), add=TRUE)
  output <- write_l1ou_placeholder_outputs(dat$tree, Y, dir=directory)
  leaf <- read.delim(output$files[["leaf"]], check.names=FALSE)
  expect_equal(leaf$node_name[leaf$param == "Y"], dat$tree$tip.label)
  expect_equal(unname(as.matrix(leaf[leaf$param == "Y", -(1:3)])), Y)
})
