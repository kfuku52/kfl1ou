small_lizard_data <- function(n_tips = 12, traits = 1, normalize = TRUE) {
  env <- new.env(parent = emptyenv())
  data("lizard.tree", package = "kfl1ou", envir = env)
  data("lizard.traits", package = "kfl1ou", envir = env)

  keep <- env$lizard.tree$tip.label[seq_len(n_tips)]
  tree <- ape::drop.tip(env$lizard.tree, setdiff(env$lizard.tree$tip.label, keep))
  Y <- as.matrix(env$lizard.traits[keep, traits, drop = FALSE])

  adjust_data(tree, Y, normalize = normalize, quietly = TRUE)
}

capture_silently <- function(expr) {
  invisible(capture.output(result <- eval.parent(substitute(expr))))
  result
}

named_input_error <- function(Y, value = 0.01) {
  stats::setNames(rep(value, nrow(Y)), rownames(Y))
}
