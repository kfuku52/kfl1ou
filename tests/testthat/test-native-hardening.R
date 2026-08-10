test_that("native pruning errors unwind cleanly and leave later calls usable", {
  effective_args <- list(
    N = 2L, n = 2L, pN = 1L, root = 3L, transa = 1,
    transb = c(0, 0), des = c(1L, 2L), anc = c(3L, 3L),
    edge = 3L
  )
  threepoint_args <- list(
    N = 2L, n = 2L, pN = 1L, dY = 1L, dX = 1L, root = 3L,
    transa = 1, transb = c(0, 0), des = c(1L, 2L),
    anc = c(3L, 3L), y = c(1, 2), X = c(1, 1)
  )

  for (index in seq_len(25L)) {
    expect_error(
      do.call(kfl1ou:::effective_sample_size_c, effective_args),
      "two or more sister"
    )
    expect_error(
      do.call(kfl1ou:::threepoint_l1ou_c, threepoint_args),
      "two or more sister"
    )
  }

  effective_args$transb <- c(1, 1)
  threepoint_args$transb <- c(1, 1)
  expect_true(all(is.finite(
    do.call(kfl1ou:::effective_sample_size_c, effective_args)
  )))
  expect_true(all(is.finite(
    do.call(kfl1ou:::threepoint_l1ou_c, threepoint_args)
  )))
})

test_that("threepoint multivariate flat workspaces agree with dense algebra", {
  set.seed(8041L)
  tree <- ape::reorder.phylo(ape::rtree(9L), "pruningwise")
  tree$root.edge <- 0.4
  response <- matrix(rnorm(18L), nrow = 9L, ncol = 2L)
  design <- matrix(rnorm(27L), nrow = 9L, ncol = 3L)

  result <- kfl1ou:::threepoint_l1ou_c(
    nrow(tree$edge), 9L, tree$Nnode, 2L, 3L, 10L,
    tree$root.edge, tree$edge.length, tree$edge[, 2L],
    tree$edge[, 1L], as.double(response), as.double(design)
  )

  covariance <- ape::vcv.phylo(tree) + tree$root.edge
  precision <- solve(covariance)
  position <- 3L
  y_one <- result[position:(position + 1L)]
  position <- position + 2L
  yy <- matrix(result[position:(position + 3L)], 2L, 2L)
  position <- position + 4L
  x_one <- result[position:(position + 2L)]
  position <- position + 3L
  xx <- matrix(result[position:(position + 8L)], 3L, 3L)
  position <- position + 9L
  xy <- matrix(result[position:(position + 5L)], 3L, 2L)

  expect_equal(
    result[[1L]],
    as.numeric(determinant(covariance, logarithm = TRUE)$modulus),
    tolerance = 1e-10
  )
  expect_equal(
    y_one,
    drop(crossprod(response, precision %*% rep(1, 9L))),
    tolerance = 1e-10
  )
  expect_equal(yy, crossprod(response, precision %*% response),
               tolerance = 1e-10)
  expect_equal(
    x_one,
    drop(crossprod(design, precision %*% rep(1, 9L))),
    tolerance = 1e-10
  )
  expect_equal(xx, crossprod(design, precision %*% design),
               tolerance = 1e-10)
  expect_equal(xy, crossprod(design, precision %*% response),
               tolerance = 1e-10)
})

test_that("group-lasso uses stable norms and rejects unrepresentable objectives", {
  scale <- 1e77
  x <- matrix(rep(scale, 2L), ncol = 1L)
  y <- rep(scale, 2L)

  lambda_max <- kfl1ou:::linreg_group_lasso_lambda_max_cpp(x, y, 1L)
  expect_true(is.finite(lambda_max))
  expect_equal(lambda_max / 4e154, 1, tolerance = 1e-12)

  path <- kfl1ou:::linreg_group_lasso_path_cpp(
    x, y, 1L, 0, tol = 1e-10
  )
  expect_true(path$converged[[1L]])
  expect_equal(path$coefficients[[1L]], 1, tolerance = 1e-10)

  too_large <- matrix(rep(1e155, 2L), ncol = 1L)
  expect_error(
    kfl1ou:::linreg_group_lasso_path_cpp(
      too_large, rep(1e155, 2L), 1L, 0
    ),
    "numerical overflow"
  )

  correlated <- matrix(1, nrow = 2L, ncol = 2L)
  stalled <- kfl1ou:::linreg_group_lasso_path_cpp(
    correlated, c(1, 1), c(1L, 1L), 0, beta_ls = 1e-31
  )
  expect_false(stalled$converged[[1L]])
  expect_equal(stalled$coefficients[, 1L], c(0, 0))
})

test_that("configuration score storage rejects NaN and remains ordered", {
  kfl1ou:::erase_configuration_score_db()
  on.exit(kfl1ou:::erase_configuration_score_db(), add = TRUE)

  expect_error(
    kfl1ou:::add_configuration_score_to_db("missing", NA_real_, ""),
    "NA or NaN"
  )
  expect_error(
    kfl1ou:::add_configuration_score_to_db("nan", NaN, ""),
    "NA or NaN"
  )

  kfl1ou:::add_configuration_score_to_db("b", 2, "second")
  kfl1ou:::add_configuration_score_to_db("a", 1, "first")
  kfl1ou:::add_configuration_score_to_db("z", Inf, "failed")
  kfl1ou:::add_configuration_score_to_db("b", 0.5, "updated")
  stored <- kfl1ou:::get_stored_config_score()

  expect_equal(stored$scores, c(0.5, 1, Inf))
  expect_equal(stored$configurations, c("b", "a", "z"))
  expect_equal(stored$moreInfo, c("updated", "first", "failed"))
  expect_equal(
    kfl1ou:::get_score_of_configuration("b"),
    list(value = 0.5, valid = TRUE)
  )
})

test_that("request-local configuration score caches are isolated", {
  first <- kfl1ou:::new_configuration_score_cache()
  second <- kfl1ou:::new_configuration_score_cache()

  kfl1ou:::add_configuration_score_to_list(c(3L, 1L), 2, "first", first)
  kfl1ou:::add_configuration_score_to_list(integer(), 1, "empty", first)
  kfl1ou:::add_configuration_score_to_list(c(1L, 3L), 4, "second", second)

  expect_equal(kfl1ou:::get_configuration_score_from_list(c(1L, 3L), first), 2)
  expect_equal(kfl1ou:::get_configuration_score_from_list(c(1L, 3L), second), 4)
  expect_true(is.na(kfl1ou:::get_configuration_score_from_list(integer(), second)))
  expect_equal(
    kfl1ou:::list_investigated_configs(first)$configurations,
    list(c(1L, 3L), integer())
  )

  kfl1ou:::erase_configuration_score_cache(first)
  expect_length(kfl1ou:::list_investigated_configs(first)$scores, 0L)
  expect_equal(kfl1ou:::get_configuration_score_from_list(c(1L, 3L), second), 4)
})

test_that("covariance kernel preserves its input edge matrix", {
  edge_list <- cbind(
    ancestor = c(3, 3, 4, 4),
    descendant = c(0, 1, 3, 2),
    length = c(1, 2, 4, 8)
  )
  original <- edge_list + 0

  result <- kfl1ou:::cmp_sqrt_OU_covariance(edge_list, 3L, 0)

  expect_identical(edge_list, original)
  expect_true(all(is.finite(result$sqrtSigma)))
  expect_true(all(is.finite(result$sqrtInvSigma)))
})

test_that("covariance kernel retains positive subnormal branch lengths", {
  tiny <- 1e-320
  skip_if(tiny == 0, "platform flushes subnormal doubles to zero")
  edge_list <- cbind(
    ancestor = c(2, 2),
    descendant = c(0, 1),
    length = c(tiny, tiny)
  )

  result <- kfl1ou:::cmp_sqrt_OU_covariance(edge_list, 2L, 0)
  reconstructed <- result$sqrtSigma %*% t(result$sqrtSigma)

  expect_true(all(is.finite(result$sqrtSigma)))
  expect_true(all(is.finite(result$sqrtInvSigma)))
  expect_gt(sum(abs(result$sqrtSigma)), 0)
  expect_equal(diag(reconstructed) / tiny, rep(1, 2L), tolerance = 1e-8)
  expect_equal(reconstructed[1L, 2L] / tiny, 0, tolerance = 1e-8)
})
