#' Adjusts the tree and traits to meet the requirements of \code{estimate_shift_configuration}
#'
#' Returns a new tree and new data matrix, where the tree edges are in
#' postorder, the data row names match the order of the tree tip labels, and
#' common pathological inputs are sanitized.
#'
#'@param tree ultrametric tree of class phylo with branch lengths.
#'@param Y trait vector/matrix.
#'@param normalize logical. If TRUE, normalizes branch lengths to a unit tree height.
#'@param quietly logical. If FALSE, changes in tree/trait are printed.
#'@param repair.tree logical. If TRUE, repairs non-positive branch lengths and
#' nearly ultrametric trees before normalization.
#'@param min.edge.length positive numeric lower bound used when repairing short
#' or non-positive branch lengths.
#'@param ultrametric.tolerance numeric vector of relative change thresholds used
#' when attempting automatic ultrametric repair.
#'@param drop.all.missing logical. If TRUE, tips with no observed trait values
#' are dropped together with their tree tips.
#'@param drop.invariant logical. If TRUE, invariant traits are removed before
#' fitting.
#'@param invariant.tolerance non-negative numeric tolerance used to detect
#' invariant traits.
#'
#'@return 
#' \item{tree}{tree of class phylo, with the same topology as the input \code{tree} but adjusted edge order.}
#' \item{Y}{trait vector/matrix with adjusted row names and row order.}
#' \item{removed.tips}{tip labels dropped because all trait values were missing.}
#' \item{removed.traits}{trait names dropped because they were invariant.}
#'@examples
#' data(lizard.tree, lizard.traits)
#' # here, lizard.traits is a matrix, so columns retain row names:
#' names(lizard.traits[,1])
#' lizard <- adjust_data(lizard.tree, lizard.traits[,1])
#' 
#' # for a data frame, make sure to retain row names if a single column is selected:
#' lizard.traits <- as.data.frame(lizard.traits)
#' lizard <- adjust_data(lizard.tree, subset(lizard.traits, select=1))
#'@export
adjust_data <- function(tree, Y, normalize = TRUE, quietly=FALSE,
                        repair.tree = TRUE,
                        min.edge.length = 1e-8,
                        ultrametric.tolerance = c(0.01, 0.05, 0.1, 1, Inf),
                        drop.all.missing = TRUE,
                        drop.invariant = TRUE,
                        invariant.tolerance = 0){

    if (!inherits(tree, "phylo"))  stop("object \"tree\" is not of class \"phylo\".")
    if (!is.null(tree$root.edge))
	    if (tree$root.edge>0)
		    stop("the tree has a non-zero root edge.")

    if (!inherits(Y, "matrix")){
        Y <- as.matrix(Y)
        if(!quietly)
            cat(paste("new Y: matrix of size", nrow(Y), "x", ncol(Y), "\n" ))
    }

    if( nrow(Y) != length(tree$tip.label)){
       stop("the number of entries/rows of the trait vector/matrix (Y) 
            doesn't match the number of tips.\n") 
    }

    if( is.null(rownames(Y)) ){
        warning("no names provided for the trait(s) entries/rows.\nAssuming that rows match the tip labels in the same order.",
                immediate.=TRUE)
        rownames(Y)  <- tree$tip.label
    } else{
        if( any(is.na(rownames(Y))) ){
            stop("some of the names in the trait vector/matrix or in the tree's 
                 tip.label are unavailable.\n")
        }
    }

    if(!identical(rownames(Y), tree$tip.label)){
        diffres = setdiff(rownames(Y), tree$tip.label)
        if( length(diffres) > 0 ){
            cat(diffres)
            stop(" do(es) not exist in the tip labels of the input tree.\n")
        }
        diffres = setdiff(tree$tip.label, rownames(Y))
        if( length(diffres) > 0 ){
            cat(diffres)
            stop(" do(es) not exist in the input trait. you may want to use
                 drop.tip(tree, setdiff(tree$tip.label,rownames(Y))) 
                 to drop extra tips in the tree.\n")
        }

        if(!quietly)
            cat("reordered the rows of the trait vector/matrix (Y) to match the order of the tip labels.\n")

        #Y  <-  Y[order(rownames(Y)),  ] 
        #Y  <-  Y[order(order(tree$tip.label)), ]
        Y <- as.matrix(Y[tree$tip.label, , drop=FALSE])

    }

    removed.tips <- character(0)
    if( drop.all.missing ){
        observed.per.tip <- rowSums(!is.na(Y))
        missing.tip.idx <- which(observed.per.tip == 0L)
        if( length(missing.tip.idx) > 0 ){
            removed.tips <- rownames(Y)[missing.tip.idx]
            if(!quietly){
                cat("dropped", length(removed.tips),
                    "tip(s) with no observed trait values.\n")
            }
            tree <- ape::drop.tip(tree, removed.tips)
            Y <- as.matrix(Y[tree$tip.label, , drop=FALSE])
        }
    }

    if( repair.tree ){
        tree <- sanitize_tree_edge_lengths(
            tree,
            min.edge.length = min.edge.length,
            quietly = quietly
        )
        tree <- repair_ultrametric_tree(
            tree,
            tolerances = ultrametric.tolerance,
            min.edge.length = min.edge.length,
            quietly = quietly
        )
    }

    if( !identical(tree$edge, reorder(tree, "postorder")$edge)){
        if(!quietly)
            cat("the new tree edges are ordered differently: in postorder.\n")
        tree  <- reorder(tree, "postorder")
    }

    if( normalize ){
        tree <- normalize_tree(tree)
        if(!quietly)
            cat("the new tree is normalized: each tip at distance 1 from the root.\n")
    }

    removed.traits <- character(0)
    if( drop.invariant ){
        filtered <- drop_invariant_traits(
            Y,
            tolerance = invariant.tolerance,
            quietly = quietly
        )
        Y <- filtered$Y
        removed.traits <- filtered$removed.traits
    }

    stopifnot(all(rownames(Y) == tree$tip.label))
    stopifnot(identical(rownames(Y), tree$tip.label))

    return(list(tree=tree, Y=Y,
                removed.tips=removed.tips,
                removed.traits=removed.traits))
}

sanitize_tree_edge_lengths <- function(tree, min.edge.length = 1e-8, quietly = FALSE){
    if( is.null(tree$edge.length) ){
        stop("the tree has no branch lengths.")
    }
    if( !is.numeric(min.edge.length) || length(min.edge.length) != 1L ||
        is.na(min.edge.length) || min.edge.length <= 0 ){
        stop("min.edge.length must be a single positive number.")
    }
    if( any(is.na(tree$edge.length)) ){
        stop("the tree has missing branch lengths.")
    }
    if( any(!is.finite(tree$edge.length)) ){
        stop("the tree has non-finite branch lengths.")
    }

    short.idx <- which(tree$edge.length < min.edge.length)
    if( length(short.idx) > 0 ){
        if(!quietly){
            cat("raised", length(short.idx),
                "branch length(s) to at least", min.edge.length, "\n")
        }
        tree$edge.length[short.idx] <- min.edge.length
    }

    tree
}

repair_ultrametric_tree <- function(tree,
                                    tolerances = c(0.01, 0.05, 0.1, 1, Inf),
                                    min.edge.length = 1e-8,
                                    quietly = FALSE){
    if( isTRUE(is.ultrametric(tree)) ){
        return(tree)
    }

    tolerances <- as.numeric(tolerances)
    tolerances <- tolerances[!is.na(tolerances) & tolerances >= 0]
    if( !length(tolerances) ){
        stop("ultrametric.tolerance must contain at least one non-negative number.")
    }

    original <- tree
    trial <- tryCatch(suppressWarnings(ape::chronoMPL(tree)), error = function(e) e)
    if( inherits(trial, "phylo") ){
        trial <- sanitize_tree_edge_lengths(
            trial,
            min.edge.length = min.edge.length,
            quietly = TRUE
        )
        sum.adjustment <- sum(abs(trial$edge.length - original$edge.length))
        total.length <- max(sum(abs(original$edge.length)), .Machine$double.eps)
        rel.change <- sum.adjustment / total.length
        accepted <- tolerances[is.infinite(tolerances) | rel.change <= tolerances]
        if( length(accepted) > 0 && isTRUE(is.ultrametric(trial)) ){
            if(!quietly){
                cat("repaired a non-ultrametric tree with chronoMPL",
                    "(relative edge-length change =", signif(rel.change, 4), ").\n")
            }
            return(trial)
        }
    }

    trial <- tryCatch(
        suppressWarnings(ape::chronos(tree, lambda = 1, quiet = TRUE)),
        error = function(e) e
    )
    if( inherits(trial, "phylo") ){
        trial <- sanitize_tree_edge_lengths(
            trial,
            min.edge.length = min.edge.length,
            quietly = TRUE
        )
        if( isTRUE(is.ultrametric(trial)) ){
            if(!quietly){
                cat("repaired a non-ultrametric tree with chronos fallback.\n")
            }
            return(trial)
        }
    }

    stop("the input tree is not ultrametric and automatic repair failed.")
}

drop_invariant_traits <- function(Y, tolerance = 0, quietly = FALSE){
    if( !is.numeric(tolerance) || length(tolerance) != 1L ||
        is.na(tolerance) || tolerance < 0 ){
        stop("invariant.tolerance must be a single non-negative number.")
    }

    n.traits <- ncol(Y)
    trait.names <- colnames(Y)
    if( is.null(trait.names) ){
        trait.names <- as.character(seq_len(n.traits))
    }

    keep <- rep(TRUE, n.traits)
    for(i in seq_len(n.traits)){
        vals <- Y[, i]
        vals <- vals[is.finite(vals)]
        if( length(vals) == 0L ){
            keep[[i]] <- FALSE
            next
        }
        keep[[i]] <- ((max(vals) - min(vals)) > tolerance)
    }

    removed.traits <- trait.names[!keep]
    if( length(removed.traits) > 0 ){
        if(!quietly){
            cat("removed", length(removed.traits),
                "invariant trait(s).\n")
        }
        Y <- Y[, keep, drop=FALSE]
    }

    if( ncol(Y) == 0L ){
        stop("all traits are invariant or entirely missing after adjustment.")
    }

    list(Y = Y, removed.traits = removed.traits)
}

restore_tip_matrix_to_tree <- function(x, original.tips, observed.tips){
    x <- as.matrix(x)
    if( is.null(rownames(x)) ){
        rownames(x) <- observed.tips
    }
    x <- x[observed.tips, , drop=FALSE]

    out <- matrix(NA_real_, nrow=length(original.tips), ncol=ncol(x))
    rownames(out) <- original.tips
    colnames(out) <- colnames(x)
    out[observed.tips, ] <- x
    out
}

edge_descendant_tip_keys <- function(tree, observed.tips = tree$tip.label){
    observed.tips <- as.character(observed.tips)
    if( !all(observed.tips %in% tree$tip.label) ){
        stop("observed.tips must be a subset of tree$tip.label.")
    }

    Z <- generate_design_matrix(tree, type="simpX")
    rownames(Z) <- tree$tip.label
    Z <- Z[observed.tips, , drop=FALSE]

    tip.sets <- lapply(seq_len(ncol(Z)), function(idx){
        sort(rownames(Z)[Z[, idx] > 0])
    })
    keys <- vapply(tip.sets, function(tips){
        if( length(tips) == 0 ){
            return("")
        }
        paste(tips, collapse="\r")
    }, character(1))

    list(tip.sets=tip.sets, keys=keys)
}

map_pruned_edges_to_original <- function(pruned.tree, original.tree,
                                         representative=c("tipward", "rootward")){
    representative <- match.arg(representative)
    pruned.tree <- reorder(pruned.tree, "postorder")
    original.tree <- reorder(original.tree, "postorder")

    if( anyDuplicated(pruned.tree$tip.label) || anyDuplicated(original.tree$tip.label) ){
        stop("tree tip labels must be unique.")
    }
    if( !all(pruned.tree$tip.label %in% original.tree$tip.label) ){
        stop("original.tree must contain all tips from the fitted tree.")
    }

    observed.tips <- pruned.tree$tip.label
    pruned.info <- edge_descendant_tip_keys(pruned.tree, observed.tips)
    original.info <- edge_descendant_tip_keys(original.tree, observed.tips)

    key.lookup <- setNames(seq_along(pruned.info$keys), pruned.info$keys)
    nonempty <- nzchar(original.info$keys)
    original.to.pruned <- rep(NA_integer_, Nedge(original.tree))
    original.to.pruned[nonempty] <- unname(key.lookup[original.info$keys[nonempty]])
    if( any(is.na(original.to.pruned[nonempty])) ){
        stop("failed to map pruned edges back to original.tree.")
    }

    pruned.to.original <- vector("list", Nedge(pruned.tree))
    for(idx in seq_len(Nedge(pruned.tree))){
        pruned.to.original[[idx]] <- which(original.to.pruned == idx)
    }

    child.depths <- tree_node_depths(original.tree)[original.tree$edge[, 2]]
    representative.edge <- vapply(pruned.to.original, function(path){
        if( length(path) == 0 ){
            return(NA_integer_)
        }
        if( representative == "tipward" ){
            return(path[[which.max(child.depths[path])]])
        }
        path[[which.min(child.depths[path])]]
    }, integer(1))

    list(
        observed.tips = observed.tips,
        original.to.pruned.edge = original.to.pruned,
        pruned.to.original.edge = pruned.to.original,
        representative.edge = representative.edge
    )
}

#' Restores a pruned fit to the original tree
#'
#' This helper is intended for reporting after a fit was run on a tree pruned
#' for missing-tip handling. It pads tip-level outputs back to the original tip
#' set and maps each shift edge to the corresponding edge path in the original
#' tree.
#'
#'@param model object of class \code{"l1ou"} returned by \pkg{kfl1ou}.
#'@param original.tree original tree before pruning tips with missing data.
#'@param representative how to choose a single representative original edge when
#' a pruned edge maps to a collapsed path in \code{original.tree}. The default
#' \code{"tipward"} picks the deepest edge on that path; \code{"rootward"}
#' picks the shallowest one.
#'
#'@return
#'An object of class \code{"restored_l1ou"} and \code{"l1ou"} with the
#' original tree, tip-level matrices padded back to the original tip set, and a
#' \code{restoration} component containing edge-path mappings.
#'
#'@export
restore_original_tree_fit <- function(model, original.tree,
                                      representative=c("tipward", "rootward")){
    if( !inherits(model, "l1ou") ){
        stop("model must inherit from class \"l1ou\".")
    }
    if( !inherits(original.tree, "phylo") ){
        stop("original.tree must be of class \"phylo\".")
    }

    representative <- match.arg(representative)
    pruned.tree <- reorder(model$tree, "postorder")
    original.tree <- reorder(original.tree, "postorder")
    mapping <- map_pruned_edges_to_original(
        pruned.tree = pruned.tree,
        original.tree = original.tree,
        representative = representative
    )

    observed.tips <- pruned.tree$tip.label
    original.tips <- original.tree$tip.label
    removed.tips <- setdiff(original.tips, observed.tips)

    restored <- model
    restored$tree <- original.tree
    restored$Y <- restore_tip_matrix_to_tree(model$Y, original.tips, observed.tips)
    restored$mu <- restore_tip_matrix_to_tree(model$mu, original.tips, observed.tips)
    restored$residuals <- restore_tip_matrix_to_tree(model$residuals, original.tips, observed.tips)
    restored$optima <- restore_tip_matrix_to_tree(model$optima, original.tips, observed.tips)

    pruned.shift.configuration <- model$shift.configuration
    shift.edge.paths <- list()
    restored.shift.configuration <- pruned.shift.configuration
    if( length(pruned.shift.configuration) > 0 ){
        pruned.idx <- as.integer(pruned.shift.configuration)
        shift.edge.paths <- mapping$pruned.to.original.edge[pruned.idx]
        names(shift.edge.paths) <- as.character(pruned.shift.configuration)
        restored.shift.configuration <- mapping$representative.edge[pruned.idx]
        names(restored.shift.configuration) <- names(pruned.shift.configuration)
    }
    restored$shift.configuration <- restored.shift.configuration

    if( !is.null(restored$l1ou.options$Z) ){
        restored$l1ou.options$Z <- generate_design_matrix(original.tree, type="simpX")
    }

    restored$restoration <- list(
        pruned.tree = pruned.tree,
        observed.tips = observed.tips,
        removed.tips = removed.tips,
        representative = representative,
        pruned.shift.configuration = pruned.shift.configuration,
        shift.edge.paths = shift.edge.paths,
        ambiguous.shifts = lengths(shift.edge.paths) > 1L,
        original.to.pruned.edge = mapping$original.to.pruned.edge,
        pruned.to.original.edge = mapping$pruned.to.original.edge
    )

    class(restored) <- c("restored_l1ou", class(model))

    if( any(restored$restoration$ambiguous.shifts) ){
        warning(
            "Some restored shifts map to collapsed edge paths on original.tree; ",
            "representative edges were chosen ", representative,
            ". Full paths are stored in $restoration$shift.edge.paths.",
            call. = FALSE
        )
    }

    restored
}

lnorm <- function(v,l=1)   { return( (sum(abs(v)^l, na.rm=TRUE))^(1/l) ) }

tree_children_edges <- function(tree){
    nNodes <- length(tree$tip.label) + tree$Nnode
    children <- vector("list", nNodes)
    for(edge.idx in seq_len(nrow(tree$edge))){
        parent <- tree$edge[edge.idx, 1]
        children[[parent]] <- c(children[[parent]], edge.idx)
    }
    return(children)
}

tree_node_depths <- function(tree, children=NULL){
    if( is.null(children) ){
        children <- tree_children_edges(tree)
    }

    nTips <- length(tree$tip.label)
    nNodes <- nTips + tree$Nnode
    root <- nTips + 1L
    depths <- rep(NA_real_, nNodes)
    depths[root] <- 0

    stack <- root
    while(length(stack) > 0){
        node <- stack[[length(stack)]]
        stack <- stack[-length(stack)]
        child.edges <- children[[node]]
        if( length(child.edges) == 0 ){
            next
        }
        for(edge.idx in child.edges){
            child <- tree$edge[edge.idx, 2]
            depths[[child]] <- depths[[node]] + tree$edge.length[[edge.idx]]
            stack <- c(stack, child)
        }
    }

    return(depths)
}

tree_edge_ages <- function(tree, children=NULL){
    if( is.null(children) ){
        children <- tree_children_edges(tree)
    }

    depths <- tree_node_depths(tree, children=children)
    nTips <- length(tree$tip.label)
    tip.height <- max(depths[seq_len(nTips)])
    tip.height - depths[tree$edge[, 1]]
}

tree_root_to_tip_edge_paths <- function(tree, children=NULL){
    if( is.null(children) ){
        children <- tree_children_edges(tree)
    }

    nTips <- length(tree$tip.label)
    root <- nTips + 1L
    paths <- vector("list", nTips)
    node.stack <- list(root)
    path.stack <- list(integer())

    while(length(node.stack) > 0){
        idx <- length(node.stack)
        node <- node.stack[[idx]]
        path <- path.stack[[idx]]
        node.stack[[idx]] <- NULL
        path.stack[[idx]] <- NULL

        child.edges <- children[[node]]
        if( length(child.edges) == 0 ){
            if( node <= nTips ){
                paths[[node]] <- path
            }
            next
        }

        for(edge.idx in rev(child.edges)){
            child <- tree$edge[edge.idx, 2]
            node.stack[[length(node.stack) + 1L]] <- child
            path.stack[[length(path.stack) + 1L]] <- c(path, edge.idx)
        }
    }

    return(paths)
}

subtree_edges_for_edge <- function(tree, edge.idx, children=NULL){
    if( is.null(children) ){
        children <- tree_children_edges(tree)
    }

    result <- integer()
    stack <- edge.idx
    while(length(stack) > 0){
        current <- stack[[length(stack)]]
        stack <- stack[-length(stack)]
        result <- c(result, current)
        child.edges <- children[[tree$edge[current, 2]]]
        if( length(child.edges) > 0 ){
            stack <- c(stack, child.edges)
        }
    }

    return(result)
}

normalize_edgelist_matrix <- function(elist){
    if( length(elist) == 0 ){
        return(matrix(character(), ncol=2))
    }
    if( is.null(dim(elist)) ){
        elist <- matrix(elist, ncol=2, byrow=TRUE)
    } else {
        elist <- as.matrix(elist)
    }
    storage.mode(elist) <- "character"
    return(elist)
}

connected_components_from_edgelist <- function(elist){
    elist <- normalize_edgelist_matrix(elist)
    if( nrow(elist) == 0 ){
        return(list())
    }

    vertices <- unique(c(elist[, 1], elist[, 2]))
    adjacency <- stats::setNames(vector("list", length(vertices)), vertices)
    for(idx in seq_len(nrow(elist))){
        u <- elist[idx, 1]
        v <- elist[idx, 2]
        adjacency[[u]] <- unique(c(adjacency[[u]], v))
        adjacency[[v]] <- unique(c(adjacency[[v]], u))
    }

    visited <- stats::setNames(rep(FALSE, length(vertices)), vertices)
    components <- list()
    for(vertex in vertices){
        if( visited[[vertex]] ){
            next
        }
        stack <- vertex
        component <- character()
        while(length(stack) > 0){
            current <- stack[[length(stack)]]
            stack <- stack[-length(stack)]
            if( visited[[current]] ){
                next
            }
            visited[[current]] <- TRUE
            component <- c(component, current)
            neighbors <- adjacency[[current]]
            if( length(neighbors) > 0 ){
                stack <- c(stack, rev(neighbors[!visited[neighbors]]))
            }
        }
        components[[length(components) + 1L]] <- component
    }

    return(components)
}

block_diag_matrix <- function(blocks){
    blocks <- Filter(function(x) !is.null(x) && length(x) > 0 && all(dim(x) > 0), blocks)
    if( length(blocks) == 0 ){
        return(matrix(0, 0, 0))
    }

    nrows <- vapply(blocks, nrow, integer(1))
    ncols <- vapply(blocks, ncol, integer(1))
    result <- matrix(0, sum(nrows), sum(ncols))

    row.offset <- 0L
    col.offset <- 0L
    for(idx in seq_along(blocks)){
        block <- as.matrix(blocks[[idx]])
        row.idx <- seq_len(nrows[[idx]]) + row.offset
        col.idx <- seq_len(ncols[[idx]]) + col.offset
        result[row.idx, col.idx] <- block
        row.offset <- row.offset + nrows[[idx]]
        col.offset <- col.offset + ncols[[idx]]
    }

    return(result)
}

l1ou_supports_multicore <- function(){
    return(.Platform$OS.type != "windows" && requireNamespace("parallel", quietly=TRUE))
}

l1ou_mclapply <- function(X, FUN, ..., mc.cores=1L){
    if( mc.cores <= 1L || !l1ou_supports_multicore() ){
        return(lapply(X, FUN, ...))
    }
    parallel_mclapply <- getFromNamespace("mclapply", "parallel")
    return(parallel_mclapply(X=X, FUN=FUN, ..., mc.cores=mc.cores))
}

l1ou_require_genlasso <- function(){
    if( !requireNamespace("genlasso", quietly=TRUE) ){
        stop("estimate_convergent_regimes(method=\"rr\") requires the optional package 'genlasso'.")
    }
}

## This function is useful for handling missing values in multivariate regression. 
## It generates a list of design matrices and trees considering according to the missing values.
gen_tree_array <- function(tree, Y){ 
    ## here I assume that the tree tip labels match the Y matrix rows 
    ## in the same order.
    tree.list <- list()
    X.1 <- generate_design_matrix(tree, type="simpX")
    rownames(X.1) <- tree$tip.label
    for(trait.idx in 1:ncol(Y)){
        availables <- rownames(Y)[!is.na(Y[,trait.idx])]

        tr <- drop.tip(tree, setdiff(tree$tip.label, availables))
        tr <- reorder(tr, "postorder")

        X.2 <- generate_design_matrix(tr, type="simpX")
        rownames(X.2) <- tr$tip.label

        ## when we drop some tips of a tree edges order changes so we need 
        ## to have a universal mapping for shift configurations.
        old.order <- rep(NA, Nedge(tree))
        for(i in 1:Nedge(tree)){
            tip.set  <- rownames(X.1)[which(X.1[,i]>0)]
            tip.set  <- intersect(tip.set, rownames(X.2)) 
            if(length(tip.set)==0)
                next
            if(length(tip.set)==1)
                edge.set <- which(X.2[tip.set,]==1)
            else
                edge.set <- which(colSums(X.2[tip.set,])==length(tip.set))

            if(length(edge.set) > 1)
                e.idx    <- edge.set[ which(colSums(X.2[,edge.set])==length(tip.set)) ]
            else
                e.idx    <- edge.set 

            stopifnot(length(e.idx)==1)
            old.order[[i]] <- e.idx
        }
        tr$old.order <- old.order
        tr$Z <- X.2
        tr$edge.age <- tree_edge_ages(tr)
        tree.list[[trait.idx]]  <-  tr
    }
    return(tree.list)
}

add_configuration_score_to_list  <- function(shift.configuration, score, moreInfo){
    shift.configuration = sort(shift.configuration)
    add_configuration_score_to_db( paste0(shift.configuration, collapse=" "), 
                                  score, moreInfo )
}

get_configuration_score_from_list <- function(shift.configuration){
    if(length(shift.configuration) > 0){
        shift.configuration <- sort(shift.configuration)
    }
    res <- get_score_of_configuration(paste0(shift.configuration, collapse=" "))
    if( res$valid == FALSE){
        return(NA)
    }
    return(res$value)
}

list_investigated_configs <- function(){
    tmpList = get_stored_config_score()
    c.s = list()
    c.s$scores = tmpList$scores
    c.s$configurations = lapply(tmpList$configurations, 
                                FUN=function(x) as.numeric(unlist(strsplit(x, split=" ")) ) ) 
    #for( i in 1:length(c.s$scores)){
    #    c.s$configurations[[i]] = as.numeric(unlist(strsplit(tmpList$configurations[[i]], split=" ")) )
    #    #c.s$moreInfo      [[i]] = as.numeric(unlist(strsplit(tmpList$moreInfo      [[i]], split=" ")) )
    #}
    return(c.s)
}

print_out <- function(eModel, quietly){
    if (quietly == FALSE){
        print( paste0( "EST: alpha: ", eModel$alpha, " sigma2: ",  
                 eModel$sigma2, " gamma: ", eModel$sigma2/(2*eModel$alpha),
                 " score: ", eModel$score, " nShifts: ", eModel$nShifts ) )
        print("-------------------")
    }
}


rescale_matrix_and_error <- function(Y, input_error=NULL){
    scale.values  <- apply(Y, 2, lnorm, l=2)
    scale.values[scale.values == 0] <- 1
    mult <- 0.1*nrow(Y)

    Y  <- mult*scale(Y, center = TRUE, scale = scale.values)
    if( !is.null(input_error) ){
        input_error <- sweep(input_error, 2, (mult/scale.values)^2, FUN="*")
    }
    return(list(Y=Y, input_error=input_error))
}

rescale_matrix <- function(Y){
    return(rescale_matrix_and_error(Y)$Y)
}



effective.sample.size <- function(phy, edges=NULL,
             model = c("BM","OUrandomRoot","OUfixedRoot","lambda","delta","EB"),
             parameters = NULL, check.pruningwise = TRUE,
             check.ultrametric = TRUE){
        # requires ultrametric tree. Kappa disallowed: causes non-ultrametric tree
        #                            both OU models result in same n_e
        # removes every edge in 'edges' to split phy into m+1 subtrees, then
        # sums log(n_e+1) over all subtrees. n_e = max(V) * one' V^{-1} one
        # where V = phylogenetic covariance matrix for the subtree.
        model = match.arg(model)
        if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
        if (check.pruningwise) phy = reorder(phy,"pruningwise")
        if (check.ultrametric)
            if (!is.ultrametric(phy))
                stop("ultrametric tree required to calculate effective sample sizes.")
        Di <- numeric(length(phy$tip.label)) # zeros
        phy <- transf.branch.lengths(phy,model,parameters=parameters,
                                     check.pruningwise=check.pruningwise,check.ultrametric=FALSE,
                                     D=Di,check.names=F)$tree
        rootedge <- dim(phy$edge)[1]+1
        if (is.null(edges)){ sortededges <- rootedge }
        else{
            o <- order(edges) 
            r <- rank(edges)
            sortededges <- c(edges[o],rootedge)
        }
        tmp <- effective_sample_size_c(
            as.integer(dim(phy$edge)[1]), # edges
            as.integer(length(phy$tip.label)), as.integer(phy$Nnode), # tips and nodes
            as.integer(length(phy$tip.label)+1), # root index
            as.double(phy$root.edge), as.double(phy$edge.length),
            as.integer(phy$edge[, 2]), as.integer(phy$edge[, 1]), # descendents and ancestors
            as.integer(sortededges) # edges to cut, including root edge
        ) # tmp has, in this order:
        if (is.null(edges))
            res <- tmp
        else res <- tmp[c(length(tmp),r)]
        return(res)
}


correct_unidentifiability <- function(tree, shift.configuration, opt){

    if( length(shift.configuration) < 2)  { return(shift.configuration) }
    shift.configuration = sort(shift.configuration)
    nTips <- nrow(opt$Z)

    all.covered.tips <- rep(FALSE, nTips)
    keep <- rep(TRUE, length(shift.configuration))
    for(idx in seq_along(shift.configuration)){
        covered.tips <- opt$Z[, shift.configuration[[idx]]] > 0
        if( !any(covered.tips & !all.covered.tips) ){
            keep[[idx]] <- FALSE
        }
        all.covered.tips <- all.covered.tips | covered.tips
    }

    if( all(keep) ){ return(shift.configuration) }

    shift.configuration <- shift.configuration[keep]
    while ( length(shift.configuration) > 1 ) {
        coverage <- rowSums(opt$Z[, shift.configuration, drop=FALSE] > 0) > 0
        if( all(coverage) ){
            shift.configuration = shift.configuration[-which.max(shift.configuration)]
        } else { break }
    }
    return(shift.configuration)
}



alpha_upper_bound <- function(tree){
    nTips       = length(tree$tip.label)
    eLenSorted  = sort(tree$edge.length[which(tree$edge[,2] < nTips)])
    eLenSorted  = eLenSorted[is.finite(eLenSorted) & (eLenSorted > 0)]
    if( !length(eLenSorted) ){
        return(Inf)
    }
    topMinLen   = min(length(eLenSorted), max(5L, ceiling(length(eLenSorted) * 0.05)))
    return( log(2)/stats::median(eLenSorted[seq_len(topMinLen)]) )
}



get_num_solutions <- function(sol.path){
    if ( grepl("lars",sol.path$call)[[1]] || grepl("l1ou_.*_path", sol.path$call)[[1]] ){
        return ( length(sol.path$beta[,1]) )
    } 

    if ( any( grepl("grplasso",sol.path$call) ) ){
        return ( length(sol.path$coefficients[1,]) )
    }
    stop(paste0(match.call(), ":undefined solver!"))
}



get_configuration_in_sol_path <- function(sol.path, index, Y, tidx=1){
    if ( grepl("lars",sol.path$call)[[1]] || grepl("l1ou_.*_path", sol.path$call)[[1]] ){
        beta      = sol.path$beta[index,]
        shift.configuration  = which( abs(beta) > 0 )
    } else if( any( grepl("grplasso",sol.path$call) ) ){
        #beta  = sol.path$coefficients[, index]
        #lIdx  = length(beta)/ncol(Y)
        #shift.configuration  = which( abs(beta[ (1+(tidx-1)*lIdx):(tidx*lIdx) ]) > 0 )
        ##NOTE: in case, in a group of variables some are zero and some non-zero i consider all as non-zero
        #shift.configuration  = which( rowSums(matrix(abs(beta),nrow=lIdx)) > 0 ) 

        beta = sol.path$coefficients[, index]
        nVariables = ncol(Y)
        MM = matrix(ifelse(abs(beta)>0,1,0), ncol = nVariables)
        shift.configuration = which(rowSums(MM) >= nVariables/2)

    } else {  
        stop(paste0(match.call(), ":undefined solver!"))
    }

    return(shift.configuration)
}



#' Converts shift values to optimum values on the edges.
#'
#' Converts a model indicated with shift values to a model with optimum values on the edges.
#'
#'@param tree ultrametric tree of class phylo with branch lengths.
#'@param shift.configuration vector of edge indices with shifts.
#'@param shift.values vector of shift values.
#'
#'@return vector of size number of edges with optimum value of the trait on the corresponding edge.
#'@examples
#' 
#' data(lizard.tree, lizard.traits)
#' 
#' sc <- c(55, 98, 118, 74, 14, 77,  32, 164)
#' sv <- c(2 ,  3,   4,  4,  1,  2, 0.5,   1)
#'
#' root.value <- -2
#'
#' optimum.values <- convert_shifts2regions(lizard.tree, sc, sv) + root.value
#'
#'
#'@export
convert_shifts2regions <-function(tree, shift.configuration, shift.values){

    stopifnot( length(shift.configuration) == length(shift.values) )

    nEdges  = Nedge(tree)
    children = tree_children_edges(tree)
    o.vec = rep(0, nEdges)

    prev.val <-options()$warn 
    options(warn = -1)
    if( length(shift.configuration) > 0)
        for(itr in 1:length(shift.configuration) ){
            eIdx     = shift.configuration[[itr]]
            o.vec.tmp = rep(0, nEdges)
            o.vec.tmp[subtree_edges_for_edge(tree, eIdx, children=children)] = shift.values[[itr]]
            o.vec = o.vec + o.vec.tmp 
        }
    options(warn = prev.val)
    return( o.vec )
}

#' Normalizes branch lengths to a unit tree height
#'
#' Normalizes all branch lengths including the root.edge if presents by the same factor, so that the distance from the root to all tips is equal to one. 
#'@param tree ultrametric tree of class phylo with branch lengths, and edges in postorder.
#'@param check.ultrametric logical. If TRUE, it checks if the input tree is ultrametric.
#'
#'@return normalized phylogenetic tree, of class phylo.
#'
#'@export
normalize_tree <- function(tree, check.ultrametric=TRUE){

    if(check.ultrametric){
        if(!is.ultrametric(tree)) 
            stop("the input tree is not ultrametric")
    }

    nTips  = length(tree$tip.label)
    children <- tree_children_edges(tree)
    depths <- tree_node_depths(tree, children=children)

    root.edge <- ifelse(is.null(tree$root.edge), 0, tree$root.edge)
    Tval     = root.edge + max(depths[seq_len(nTips)])
    #Tval = mean ( sapply( 1:nTips, FUN=function(x) sum(tree$edge.length[root2tip[[x]]])   )  )

    tree$edge.length = tree$edge.length / Tval
    if(!is.null(tree$root.edge)){
	    tree$root.edge <- tree$root.edge / Tval 
    }
    return(tree)
}


#'
#' Visualizes a shift configuration: tree and trait(s)
#'
#' Plots the tree annotated to show the edges with a shift, and the associated trait data side by side.
#' Each trait is shown standardized, to visually highlight the species with low,
#' average or high values. In other words, the axis scale shows the values of
#' (trait - m)/sd where m is the observed trait mean and sd is the observed trait
#' standard deviation. The values on the left and right side of the axis are the
#' mininum and maximum of the standardized trait: (min-m)/sd and (max-m)/sd.
#' The place where the bars start corresponds to the mean of the original trait values.
#'
#'@param x object of class l1ou returned by \pkg{kfl1ou}.
#'@param palette vector of colors, of size the number of shifts plus one. The last element is the color for the background regime (regime at the root).
#'@param edge.shift.ann logical. If TRUE, annotates edges by shift values. 
#'@param edge.shift.adj adjustment argument to give to edgelabel() for labeling edges by shift values.
#'@param asterisk logical. If TRUE, the shift positions will be annotated by "*". It is useful for gray scale plots.
#'@param edge.label vector of size number of edges.
#'@param edge.label.ann logical. If TRUE, annotates edges by labels in tree$edge.label, if non-empty, or edge.label. 
#'@param edge.label.adj adjustment argument to give to edgelabel() for labeling edges.
#'@param edge.label.pos relative position of the edge.label on the edge. 0 for the beginning of the edge and 1 for the end of the edge. 
#'@param edge.ann.cex amount by which the annotation text should be magnified relative to the default.
#'@param plot.bar logical. If TRUE, the bars corresponding to the trait values will be plotted.
#'@param bar.axis logical. If TRUE, the axis of of trait(s) range will be plotted. 
#'@param ... further arguments to be passed on to plot.phylo. 
#'
#'@return none.
#'@examples
#' 
#' data(lizard.traits, lizard.tree)
#' keep <- lizard.tree$tip.label[1:15]
#' tree <- drop.tip(lizard.tree, setdiff(lizard.tree$tip.label, keep))
#' tree <- reorder(tree, "postorder")
#' Y <- lizard.traits[keep, 1]
#' eModel <- estimate_shift_configuration(tree, Y, criterion="AICc", max.nShifts=2)
#' nEdges <- Nedge(tree)
#' ew <- rep(1,nEdges) 
#' if (length(eModel$shift.configuration) > 0) {
#'   ew[eModel$shift.configuration] <- 3
#' }
#' plot(eModel, cex=0.5, label.offset=0.02, edge.width=ew)
#'
#'@export
#'
plot.l1ou <- function (x, palette = NA, 
                       edge.shift.ann=TRUE,  edge.shift.adj=c(0.5,-.025),
                       edge.label=c(), asterisk = TRUE,
                       edge.label.ann=FALSE, edge.label.adj=c(0.5,    1), 
                       edge.label.pos=NA,
                       edge.ann.cex = 1, 
                       plot.bar = TRUE, bar.axis = TRUE, ...) 
{
    model = x
    tree = model$tree
    s.c = model$shift.configuration
    stopifnot(identical(tree$edge, reorder(tree, "postorder")$edge))
    nShifts = model$nShifts
    nEdges = Nedge(tree)
    if (bar.axis) 
        par(oma = c(3, 0, 0, 3))
    Y = as.matrix(model$Y)
    stopifnot(identical(rownames(Y), tree$tip.label))

    if (bar.axis) 
        par(oma = c(3, 0, 0, 3))

    if (plot.bar) {
        layout(matrix(c(1+ncol(Y),1:ncol(Y)), nrow=1), 
               widths = c(2,rep(1, ncol(Y)))
               )
    }

    #NOTE: assiging colors the shifts
    if (all(is.na(palette))) {
        palette = c(sample(rainbow(nShifts)), "gray")
        if( !is.null(names(s.c)) ){
            ids = unique(names(s.c))
            tmp = sample(rainbow(length(ids)))
	    for( id in ids ){
		    if(id == "0"){
			    palette[which(names(s.c) == id)] = "gray"  #background 
			    next
		    }
		    palette[which(names(s.c)==id)] = tmp[which(ids==id)]
	    }
        }
    }

    stopifnot(length(palette) == model$nShifts + 1)

    edgecol = rep(palette[nShifts + 1], nEdges)
    counter = 1
    Z = model$l1ou.options$Z
    if(length(s.c) > 0)
        for (shift in sort(s.c, decreasing = T)) {
            edgecol[[shift]] = palette[[which(s.c == shift)]]
            tips = which(Z[, shift] > 0)
            for (tip in tips) {
                edgecol[which(Z[tip, 1:shift] > 0)] = palette[[which(s.c == 
                    shift)]]
            }
            counter = counter + 1
        }



    ##A dummy plot just to get the plotting order
    plot.phylo(tree, plot=FALSE)
    ape_plot_env = getFromNamespace(".PlotPhyloEnv", "ape")
    lastPP = get("last_plot.phylo", envir = ape_plot_env)
    o = order(lastPP$yy[1:length(tree$tip.label)])
    par.new.default <- par()$new ##just to be careful with the global variable
    par(new=TRUE)

    #NOTE: plotting bar plot .....
    if (plot.bar) {
        nTips = length(tree$tip.label)
        barcol = rep("gray", nTips)
        for (i in 1:nTips) {
            barcol[[i]] = edgecol[which(tree$edge[, 2] == i)]
        }
        if (bar.axis) 
            par(mar = c(0, 0, 0, 3))
        for (i in 1:ncol(Y)) {
            normy = (Y[, i] - mean(Y[, i], na.rm=TRUE))/sd(Y[, i], na.rm=TRUE)
            barplot(as.vector(normy[o]), border = FALSE, col = barcol[o], 
                    horiz = TRUE, names.arg = "", xaxt = "n")
            if (bar.axis){
                axis(1, at = range(normy, na.rm=TRUE), 
                     labels = round(range(normy, na.rm=TRUE), 
                                    digits = 2))
            }
            if (!is.null(colnames(Y)) && length(colnames(Y)) > 
                (i - 1)) 
                mtext(colnames(Y)[[i]], cex = 1, line = +1, side = 1)
        }
    }

    #NOTE: plotting the tree etc etc
    plot.phylo(tree, edge.color = edgecol, no.margin = TRUE, 
        ...)

    if (length(s.c) > 0) {
        if (asterisk) {
            Z = generate_design_matrix(tree, type = "apprX")
            for (idx in 1:length(s.c)) {
                sP = s.c[[idx]]
                pos = max(Z[, sP])
                edge.labels = rep(NA, length(tree$edge[, 1]))
                edge.labels[sP] = "*"
                edgelabels(edge.labels, cex = 3 * edge.ann.cex, 
                  adj = c(0.5, 0.8), frame = "none", date = pos)
            }
        }
    }
    if (edge.shift.ann) {
        eLabels = rep(NA, nEdges)
        for (shift in s.c) {
            eLabels[shift] = paste(round(model$shift.values[which(s.c == 
                shift), ], digits = 2), collapse = ",")
        }
        edgelabels(eLabels, cex = edge.ann.cex, adj = edge.shift.adj, 
            frame = "none")
    }
    if (edge.label.ann) {
        if (length(tree$edge.label) == 0) {
            if (length(edge.label) == 0) {
                stop("no edge labels are provided via tree$edge.label or edge.label!")
            }
            tree$edge.label = edge.label
        }
        Z = generate_design_matrix(tree, type = "apprX")
        if (!is.na(edge.label.pos)) 
            if (edge.label.pos < 0 || edge.label.pos > 1) 
                stop("edge.label.pos should be between 0 and 1")
        for (idx in 1:length(tree$edge.label)) {
            if (is.na(tree$edge.label[[idx]])) 
                next
            pos = max(Z[, idx])
            if (!is.na(edge.label.pos)) {
                pos = pos - edge.label.pos * tree$edge.length[[idx]]
            }
            edge.labels = rep(NA, length(tree$edge[, 1]))
            edge.labels[[idx]] = tree$edge.label[[idx]]
            edgelabels(edge.labels, cex = edge.ann.cex, adj = edge.label.adj, 
                frame = "none", date = pos)
        }
    }
    par(new=par.new.default)
}

#'@export
plot.restored_l1ou <- function(x, ...){
    tmp <- x
    class(tmp) <- setdiff(class(tmp), "restored_l1ou")
    plot.l1ou(tmp, ...)
}

#'
#' Prints out a summary of the shift configurations investigated by \code{\link{estimate_shift_configuration}}  
#'
#' prints the list of the shift configurations sorted by number of shifts and corresponding ic scores.
#'
#'@param fitted object of class l1ou returned by \pkg{kfl1ou}.
#'@param ... further arguments. 
#'
#'@return 
#'\item{shift.configurations}{list of shift configurations sorted by number of shifts.}
#'\item{scores}{numeric vector of scores corresponding to shift.configurations.}
#'\item{nShifts}{integer vector giving the number of shifts in each configuration.}
#'
#'@examples
#' 
#' data(lizard.traits, lizard.tree)
#' keep <- lizard.tree$tip.label[1:15]
#' tree <- drop.tip(lizard.tree, setdiff(lizard.tree$tip.label, keep))
#' tree <- reorder(tree, "postorder")
#' Y <- lizard.traits[keep, 1]
#' eModel <- estimate_shift_configuration(tree, Y, criterion="AICc", max.nShifts=2)
#' model.profile  <- profile(eModel)
#' plot(model.profile$nShifts, model.profile$scores)
#'
#'@export
#'
profile.l1ou <- function(fitted, ...)
{
    model <- fitted

    profile.data = model$profile
    if( is.null(profile.data) || is.null(profile.data$scores) ||
        is.null(profile.data$configurations) || length(profile.data$scores) == 0L ){
        if( is.null(model$score) || is.null(model$shift.configuration) ){
            stop("profile information is unavailable for this model.")
        }
        profile.data <- list(scores = model$score,
                             configurations = list(model$shift.configuration))
    }

    if( length(profile.data$scores) != length(profile.data$configurations) ){
        stop("model$profile is malformed: scores and configurations must have the same length.")
    }

    p.d = list(shift.configurations = list(), scores = numeric(), nShifts = integer())
    lens = vapply(profile.data$configurations, length, integer(1))
    ord = order(lens, profile.data$scores)
    profile.data$scores = profile.data$scores[ord]
    profile.data$configurations = profile.data$configurations[ord]

    clength = NA_integer_
    counter = 0L
    for (i in seq_along(profile.data$scores)) {
        if (!is.na(clength) && clength == length(profile.data$configurations[[i]])) {
            next
        }
        clength = length(profile.data$configurations[[i]])
        counter = counter + 1L
        p.d$shift.configurations[[counter]] = profile.data$configurations[[i]]

        p.d$nShifts[counter] = length(profile.data$configurations[[i]])
        p.d$scores[counter] = profile.data$scores[[i]]
        #p.d$gamma  [[counter]] = profile.data$moreInfo[[i]][[1]] ##the stationary variance
        #p.d$logLik [[counter]] = profile.data$moreInfo[[i]][[2]] ##the log likelihood
    }
    return(p.d)
}

#'
#' Returns the best shift configuration with a given number of shifts among the shift configurations that have been evaluated.
#'
#'@param model object of class l1ou returned by \pkg{kfl1ou}.
#'@param nShifts number of shifts.
#'
#'@return indices of the edges with shifts
#'
#'@export
get_shift_configuration <- function(model, nShifts){
    p.d = profile(model) 
    if( nShifts > length(p.d$shift.configurations)+1) # starts at 0 shifts
        stop("There is no configuration with the given number of shifts")

    for( i in seq_along(p.d$shift.configurations)){
        if( length(p.d$shift.configurations[[i]]) == nShifts)
            return(p.d$shift.configurations[[i]])
    }
    stop("There is no configuration with the given number of shifts")
}

#'
#' Prints out a summary of the model 
#'
#' prints out a summary of the model 
#'
#'@param model object of class l1ou returned by \pkg{kfl1ou}.
#'@param nTop.scores number of top scores and shift configuration to print out.
#'@param ... further arguments. 
#'
#'@return none.
#'@examples
#' 
#' data(lizard.traits, lizard.tree)
#' keep <- lizard.tree$tip.label[1:15]
#' tree <- drop.tip(lizard.tree, setdiff(lizard.tree$tip.label, keep))
#' tree <- reorder(tree, "postorder")
#' Y <- lizard.traits[keep, 1]
#' eModel <- estimate_shift_configuration(tree, Y, criterion="AICc", max.nShifts=2)
#' summary(eModel)
#'
#'@export
#'
summary.l1ou <- function(object, nTop.scores=5, ...){
    model <- object
    cat("number of shifts: ")
    cat(model$nShifts)
    cat("\n")

    cat("edge indices of the shift configuration (column names) and the corresponding shift values:\n")
    tmp.mat = t(as.matrix(model$shift.values))
    if(length(model$shift.configuration)>0)
        colnames(tmp.mat) = model$shift.configuration
    if(!all(is.null(colnames(model$Y)))){
        rownames(tmp.mat) = colnames(model$Y)
    }
    print(tmp.mat)

    cat("shift edges and corresponding jump in the trait means at the tips:")
    if(length(model$shift.configuration) > 0){
      tmp.mat = t(as.matrix(model$shift.means))
      colnames(tmp.mat) = model$shift.configuration
      if(!all(is.null(colnames(model$Y)))){
        rownames(tmp.mat) = colnames(model$Y)
      }
      print(tmp.mat)
    } else {
      cat("\n")
    }

    cat("\n")
    cat(paste0(model$l1ou.options$criterion, " score: "))
    cat(model$score)
    cat("\n")

    tmp.mat = rbind(model$alpha, 
                    model$sigma2, 
                    model$sigma2/(2 * model$alpha))
    rownames(tmp.mat) = c("adaptation rate (alpha)", 
                          "variance (sigma2)", 
                          "stationary variance (gamma)")
    if( !is.null(model$sigma2_error) &&
        any(is.finite(model$sigma2_error) & model$sigma2_error > 0) ){
        tmp.mat = rbind(tmp.mat, model$sigma2_error)
        rownames(tmp.mat)[nrow(tmp.mat)] = "measurement error variance (sigma2_error)"
    }
    tmp.mat = rbind(tmp.mat, model$logLik)
    rownames(tmp.mat)[nrow(tmp.mat)] = "logLik"
    if(!all(is.null(colnames(model$Y)))){
        colnames(tmp.mat) = colnames(model$Y)
    }
    print(tmp.mat)
    cat("\n")

    #cat("\n")
    #cat("optimum values at tips: \n")
    #print(model$optima) # too long vector: # of tips, could be thousands
    #cat("\nexpected values at the tips:\n")
    #print(model$mu)

    top.scores = min(nTop.scores, length(model$profile$scores))
    if (top.scores>0){
    cat(paste0(c("\ntop", top.scores, "best scores among candidate models evaluated during the search:\n")))
    cat("scores\t\tshift.configurations\n")
    for (i in 1:top.scores){
        cat(model$profile$scores[[i]])
        cat("\t")
        cat(model$profile$configurations[[i]])
        cat("\n")
    }
    }
}

#'@export
print.l1ou <- function(x, ...){
    model <- x
    cat("number of shifts: ")
    cat(model$nShifts)
    cat("\n")

    cat(paste0(model$l1ou.options$criterion, " score: "))
    cat(model$score)
    if(!is.null(model$cr.score)){
        cat("\n")
        cat(paste0(model$l1ou.options$criterion, " CR score: "))
        cat(model$cr.score)
    }
    cat("\n")

    cat("edge indices of the shift configuration (column names) and the corresponding shift values:\n")

    if( length(model$shift.configuration) > 0){
        tmp.mat = t(as.matrix(model$shift.values))
        if(length(model$shift.configuration)>0)
            colnames(tmp.mat) = model$shift.configuration
        if(!all(is.null(colnames(model$Y)))){
            rownames(tmp.mat) = colnames(model$Y)
        }
        print(tmp.mat)
        cat("\n")
    }


    sc <- model$shift.configuration
    if( !is.null(names(sc)) ){
        cat("convergent regimes and edge indices of the shift configuration\n")
        for( reg in sort(unique(names(sc))) )
            cat( paste0( "regime ", reg, "-> ", paste0(sc[which(names(sc)==reg)], collapse=", "), "\n" ) )
        cat("\n")
    }


    tmp.mat <- rbind(model$alpha, 
                    model$sigma2, 
                    model$sigma2/(2 * model$alpha))
    rownames(tmp.mat) <- c("adaptation rate (alpha)", 
                          "variance (sigma2)", 
                          "stationary variance (gamma)")
    if( !is.null(model$sigma2_error) &&
        any(is.finite(model$sigma2_error) & model$sigma2_error > 0) ){
        tmp.mat <- rbind(tmp.mat, model$sigma2_error)
        rownames(tmp.mat)[nrow(tmp.mat)] <- "measurement error variance (sigma2_error)"
    }
    tmp.mat <- rbind(tmp.mat, model$logLik)
    rownames(tmp.mat)[nrow(tmp.mat)] <- "logLik"
    if(!all(is.null(colnames(model$Y)))){
        colnames(tmp.mat) <- colnames(model$Y)
    }
    print(tmp.mat)
    cat("\n")
}

#'@export
print.restored_l1ou <- function(x, ...){
    tmp <- x
    class(tmp) <- setdiff(class(tmp), "restored_l1ou")
    print.l1ou(tmp, ...)

    removed.tips <- x$restoration$removed.tips
    if( length(removed.tips) > 0 ){
        cat("restoration note: padded", length(removed.tips),
            "dropped tip(s) back onto the original tree.\n")
    }
    if( any(x$restoration$ambiguous.shifts) ){
        cat("restoration note: some shifts remain ambiguous on the original tree; ",
            "representative edges were chosen ", x$restoration$representative,
            ". Full paths are stored in $restoration$shift.edge.paths.\n",
            sep="")
    }
    invisible(x)
}
