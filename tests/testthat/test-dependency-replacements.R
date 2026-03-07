test_that("block_diag_matrix places blocks on the diagonal", {
  block1 <- matrix(1:4, nrow = 2)
  block2 <- matrix(5:10, nrow = 3)

  out <- kfl1ou:::block_diag_matrix(list(block1, block2))

  expect_equal(dim(out), c(5, 4))
  expect_equal(out[1:2, 1:2], block1)
  expect_equal(out[3:5, 3:4], block2)
  expect_equal(out[1:2, 3:4], matrix(0, 2, 2))
  expect_equal(out[3:5, 1:2], matrix(0, 3, 2))
})

test_that("connected_components_from_edgelist finds undirected components", {
  elist <- rbind(
    c("0", "0"),
    c("1", "2"),
    c("2", "3"),
    c("4", "4")
  )

  components <- kfl1ou:::connected_components_from_edgelist(elist)
  components <- lapply(components, sort)

  expect_equal(components[[1]], "0")
  expect_equal(components[[2]], c("1", "2", "3"))
  expect_equal(components[[3]], "4")
})

test_that("generate_design_matrix matches explicit root-to-tip paths without igraph", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  tree <- ape::reorder.phylo(tree, "postorder")

  X <- kfl1ou:::generate_design_matrix(tree, type = "simpX")
  rownames(X) <- tree$tip.label

  parent_edge <- integer(length(tree$tip.label) + tree$Nnode)
  parent_node <- integer(length(tree$tip.label) + tree$Nnode)
  for (edge_idx in seq_len(nrow(tree$edge))) {
    parent_node[tree$edge[edge_idx, 2]] <- tree$edge[edge_idx, 1]
    parent_edge[tree$edge[edge_idx, 2]] <- edge_idx
  }

  for (tip_idx in seq_along(tree$tip.label)) {
    expected <- integer()
    node <- tip_idx
    while (node != length(tree$tip.label) + 1L) {
      expected <- c(parent_edge[node], expected)
      node <- parent_node[node]
    }
    expect_equal(sort(which(X[tip_idx, ] == 1)), sort(expected))
  }
})

test_that("convert_shifts2regions propagates shifts through descendant edges", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  tree <- ape::reorder.phylo(tree, "postorder")
  X <- kfl1ou:::generate_design_matrix(tree, type = "simpX")
  rownames(X) <- tree$tip.label

  left_clade <- which(
    colSums(X[c("a", "b"), , drop = FALSE]) == 2 &
      colSums(X[c("c", "d"), , drop = FALSE]) == 0
  )
  tip_a <- which(X["a", ] == 1 & colSums(X) == 1)

  expect_length(left_clade, 1)
  expect_length(tip_a, 1)

  regions <- convert_shifts2regions(tree, c(left_clade, tip_a), c(2, -1))

  expect_equal(regions[left_clade], 2)
  expect_equal(regions[tip_a], 1)
  expect_equal(regions[which(X["b", ] == 1 & colSums(X) == 1)], 2)
  expect_equal(regions[which(X["c", ] == 1 & colSums(X) == 1)], 0)
  expect_equal(regions[which(X["d", ] == 1 & colSums(X) == 1)], 0)
})
