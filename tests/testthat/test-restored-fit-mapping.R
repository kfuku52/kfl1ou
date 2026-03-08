test_that("restore_original_tree_fit pads dropped tips back onto the original tree", {
  env <- new.env(parent = emptyenv())
  data("lizard.tree", package = "kfl1ou", envir = env)
  data("lizard.traits", package = "kfl1ou", envir = env)

  keep <- env$lizard.tree$tip.label[seq_len(8)]
  original.tree <- ape::drop.tip(env$lizard.tree, setdiff(env$lizard.tree$tip.label, keep))
  Y <- as.matrix(env$lizard.traits[keep, 1, drop = FALSE])
  rownames(Y) <- keep
  Y[1, 1] <- NA_real_

  dat <- capture_silently(adjust_data(original.tree, Y, normalize = FALSE, quietly = FALSE))
  fit <- suppressWarnings(
    capture_silently(
      estimate_shift_configuration(
        dat$tree,
        dat$Y,
        max.nShifts = 1,
        criterion = "AICc",
        quietly = TRUE
      )
    )
  )
  restored <- restore_original_tree_fit(fit, original.tree)

  expect_s3_class(restored, "restored_l1ou")
  expect_s3_class(restored, "l1ou")
  expect_equal(restored$restoration$removed.tips, keep[[1]])
  expect_true(all(is.na(restored$Y[keep[[1]], ])))
  expect_true(all(is.na(restored$mu[keep[[1]], ])))
  expect_equal(restored$Y[fit$tree$tip.label, , drop = FALSE], fit$Y)
  expect_equal(
    unname(restored$mu[fit$tree$tip.label, , drop = FALSE]),
    unname(fit$mu)
  )
})

test_that("restore_original_tree_fit records collapsed edge paths", {
  original.tree.raw <- ape::read.tree(text = "(((A:1,B:1):1,C:2):1,D:3);")
  Y <- matrix(c(3, NA, 0, 0), ncol = 1, dimnames = list(original.tree.raw$tip.label, "trait1"))

  dat <- capture_silently(adjust_data(original.tree.raw, Y, normalize = FALSE, quietly = FALSE))
  tip.edge.pruned <- which(dat$tree$edge[, 2] == match("A", dat$tree$tip.label))
  fit <- suppressWarnings(
    capture_silently(
      fit_OU(
        dat$tree,
        dat$Y,
        shift.configuration = tip.edge.pruned,
        criterion = "AICc"
      )
    )
  )

  original.tree <- stats::reorder(original.tree.raw, "postorder")
  restored <- suppressWarnings(restore_original_tree_fit(fit, original.tree.raw, representative = "tipward"))
  shift.path <- restored$restoration$shift.edge.paths[[1]]
  original.tip.edge <- which(original.tree$edge[, 2] == match("A", original.tree$tip.label))
  removed.tip.edge <- which(original.tree$edge[, 2] == match("B", original.tree$tip.label))

  expect_length(shift.path, 2)
  expect_true(original.tip.edge %in% shift.path)
  expect_equal(restored$shift.configuration[[1]], original.tip.edge)
  expect_true(is.na(restored$restoration$original.to.pruned.edge[removed.tip.edge]))
})
