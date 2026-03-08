test_that("write_l1ou_placeholder_outputs writes expected placeholder files", {
  dat <- small_lizard_data(n_tips = 8, traits = 1:2)
  Y <- dat$Y
  Y[] <- 1

  out.dir <- file.path(tempdir(), paste0("kfl1ou-placeholder-", Sys.getpid()))
  unlink(out.dir, recursive = TRUE, force = TRUE)

  out <- write_l1ou_placeholder_outputs(dat$tree, Y, dir = out.dir, reason = "all invariant")

  expect_true(all(file.exists(unname(out$files))))
  expect_equal(out$reason, "all invariant")

  tree.tbl <- read.delim(out$files[["tree"]], check.names = FALSE)
  regime.tbl <- read.delim(out$files[["regime"]], check.names = FALSE)
  leaf.tbl <- read.delim(out$files[["leaf"]], check.names = FALSE)

  expect_equal(names(tree.tbl),
               c("num_shift", "num_regime", "num_conv_regime",
                 "num_uniq_regime", "num_species", "num_leaf", "model_score"))
  expect_equal(tree.tbl$num_shift, 0)
  expect_equal(regime.tbl$param, c("alpha", "sigma2", "intercept", "log_likelihood"))
  expect_equal(sort(unique(leaf.tbl$param)), c("Y", "mu", "optima", "residuals"))
  expect_equal(sort(unique(leaf.tbl$node_name)), sort(dat$tree$tip.label))
  expect_true(file.info(out$files[["plot"]])$size > 0)
})
