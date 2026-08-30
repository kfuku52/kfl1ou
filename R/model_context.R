# Keep reusable settings separate from data-dependent working state. Older
# fitted objects can contain these caches; discard them at a public boundary.
l1ou_model_options <- function(opt){
    cached <- c(
        "Z", "design.edges", "edge.age", ".tree.identity", "tree.list",
        "prepared.tree", "prepared.tree.list", "multivariate.tree.cache",
        "multivariate.pruning.cache", "score.cache", "search.state",
        "effective.search.strategy", "sparse.reference.configuration",
        "ageMatrix", "shift.configuration"
    )
    opt[cached] <- NULL
    opt
}

l1ou_tree_identity <- function(tree){
    list(edge=tree$edge, edge.length=tree$edge.length,
         tip.label=tree$tip.label, Nnode=tree$Nnode,
         scale=l1ou_tree_time_scale(tree))
}

l1ou_design_columns <- function(tree, opt, edges){
    edges <- as.integer(edges)
    if(!is.null(opt$Z)){
        columns <- if(is.null(opt$design.edges)) edges else
            match(edges, opt$design.edges)
        if(!anyNA(columns)) return(opt$Z[, columns, drop=FALSE])
    }
    generate_design_matrix(tree, "simpX", edges=edges)
}

# Apply explicit overrides to the original settings, including options for
# candidate generation. Full selection repeats convergence discovery when the
# input model was convergent; conditional selection preserves its partition.
refit_l1ou_model <- function(model, Y=model$Y,
                             selection=c("conditional", "full"),
                             overrides=list()){
    selection <- match.arg(selection)
    opt <- make_nested_l1ou_options(l1ou_model_options(model$l1ou.options))
    if(selection == "conditional" &&
       isTRUE(opt$convergent.fixed.alpha) &&
       !any(c("alpha.lower.bound", "alpha.upper.bound") %in% names(overrides))){
        opt$alpha.lower.bound <- as.numeric(model$alpha)
        opt$alpha.upper.bound <- as.numeric(model$alpha)
        opt$alpha.starting.value <- as.numeric(model$alpha)
    }
    opt[names(overrides)] <- overrides
    opt$use.saved.scores <- identical(selection, "full")
    if(selection == "conditional"){
        return(fit_OU(
            model$tree, Y, as.integer(model$shift.configuration),
            cr.regimes=model$convergent.regimes, l1ou.options=opt
        ))
    }
    fit <- estimate_shift_configuration(model$tree, Y, l1ou.options=opt)
    if(isTRUE(model$convergent)){
        method <- model$l1ou.options$method
        if(is.null(method)) method <- "backward"
        fit <- estimate_convergent_regimes(
            fit, criterion=opt$criterion, method=method,
            fixed.alpha=isTRUE(model$l1ou.options$convergent.fixed.alpha),
            nCores=1L
        )
    }
    fit
}

l1ou_model_mean_design_blocks <- function(model){
    Y <- as.matrix(model$Y)
    shifts <- as.integer(model$shift.configuration)
    alpha <- rep(as.numeric(model$alpha), length.out=ncol(Y))
    if(!isTRUE(model$convergent) || is.null(model$convergent.regimes)){
        base <- cbind(intercept=1, generate_design_matrix(
            model$tree, "simpX", edges=shifts
        ))
        return(rep(list(base), ncol(Y)))
    }
    # Marginal fits with missing data use pruned trait trees and mapped shifts.
    trait.trees <- if(anyNA(Y) && is.null(model$trait.covariance)){
        gen_tree_array(model$tree, Y, edges=shifts)
    } else NULL
    lapply(seq_len(ncol(Y)), function(j){
        tree <- model$tree
        trait.shifts <- shifts
        regimes <- model$convergent.regimes
        if(!is.null(trait.trees)){
            tree <- trait.trees[[j]]
            mapped <- tree$old.order[shifts]
            trait.shifts <- as.integer(mapped[!is.na(mapped)])
            regimes <- map_convergent_regimes_to_trait(regimes, shifts, mapped)
        }
        design <- generate_prediction_vec(
            tree, trait.shifts, regimes, alpha[[j]],
            root.model=model$l1ou.options$root.model
        )
        result <- matrix(0, nrow(Y), ncol(design))
        result[match(tree$tip.label, model$tree$tip.label), ] <- design
        result
    })
}
