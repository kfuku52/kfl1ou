#' Tip partition induced by a shift configuration
#'
#' Converts exact edge placements into the partition of extant tips induced by
#' their ordered ancestry through the selected shift edges. The integer labels
#' are canonical and therefore do not depend on arbitrary regime names.
#'
#'@param tree phylogenetic tree.
#'@param shift.configuration integer shift-edge indices.
#'@return A named integer vector with one canonical regime label per tip.
#'@export
shift_tip_partition <- function(tree, shift.configuration=integer()){
    validate_l1ou_tree(tree, require.positive.edges=FALSE)
    shift.configuration <- sort(unique(as.integer(shift.configuration)))
    if(anyNA(shift.configuration) ||
       any(shift.configuration < 1L | shift.configuration > Nedge(tree))){
        stop("shift.configuration contains an invalid edge index.")
    }
    if(length(shift.configuration) == 0L){
        return(stats::setNames(rep(1L, length(tree$tip.label)), tree$tip.label))
    }
    Z <- generate_design_matrix(tree, "simpX")[
        , shift.configuration, drop=FALSE
    ] > 0
    ancestry <- apply(Z, 1L, function(row) paste(as.integer(row), collapse=""))
    clusters <- split(tree$tip.label, ancestry)
    cluster.keys <- vapply(clusters, function(tips){
        paste(sort(tips), collapse="\r")
    }, character(1))
    ordered.keys <- names(sort(cluster.keys))
    labels <- match(ancestry, ordered.keys)
    stats::setNames(as.integer(labels), tree$tip.label)
}


#' Canonical key for a shift-induced tip partition
#'
#'@inheritParams shift_tip_partition
#'@return A canonical character representation of the unordered tip groups.
#'@export
shift_partition_key <- function(tree, shift.configuration=integer()){
    partition <- shift_tip_partition(tree, shift.configuration)
    groups <- split(names(partition), partition)
    keys <- sort(vapply(groups, function(tips){
        paste(sort(tips), collapse="\r")
    }, character(1)))
    paste(keys, collapse="\n")
}


partition_coassignment_matrix <- function(partition){
    partition <- as.integer(partition)
    outer(partition, partition, `==`) + 0
}


adjusted_rand_index <- function(left, right){
    left <- as.integer(left)
    right <- as.integer(right)
    if(length(left) != length(right)) stop("partitions must have equal lengths.")
    n <- length(left)
    if(n < 2L) return(1)
    table.lr <- table(left, right)
    choose2 <- function(x) x * pmax(x - 1, 0) / 2
    index <- sum(choose2(table.lr))
    left.sum <- sum(choose2(rowSums(table.lr)))
    right.sum <- sum(choose2(colSums(table.lr)))
    total <- choose2(n)
    expected <- left.sum * right.sum / total
    maximum <- 0.5 * (left.sum + right.sum)
    denominator <- maximum - expected
    if(abs(denominator) <= .Machine$double.eps){
        return(if(identical(as.integer(left), as.integer(right))) 1 else 0)
    }
    (index - expected) / denominator
}


#' Enumerate shift configurations equivalent at the tips
#'
#' Enumerates a bounded configuration space and returns all parsimonious
#' configurations that induce the same unordered partition of extant tips.
#'
#'@inheritParams shift_tip_partition
#'@param max.nShifts maximum number of shifts to enumerate.
#'@param candid.edges optional allowed edge indices.
#'@param max.configurations enumeration safety limit.
#'@return A list of equivalent integer shift configurations.
#'@export
equivalent_shift_configurations <- function(
        tree, shift.configuration, max.nShifts=length(shift.configuration),
        candid.edges=NA, max.configurations=5000L){
    validate_l1ou_tree(tree, require.positive.edges=FALSE)
    opt <- normalize_shift_search_options(
        list(
            max.nShifts=max.nShifts,
            candid.edges=candid.edges,
            exhaustive.max.configurations=max.configurations,
            search.strategy="exhaustive",
            Z=generate_design_matrix(tree, "simpX")
        ),
        tree,
        matrix(0, nrow=length(tree$tip.label), ncol=1L,
               dimnames=list(tree$tip.label, "trait"))
    )
    resolve_shift_search_strategy(tree, opt)
    target <- shift_partition_key(tree, shift.configuration)
    configurations <- enumerate_shift_configurations(tree, opt)
    configurations[vapply(configurations, function(configuration){
        identical(shift_partition_key(tree, configuration), target)
    }, logical(1))]
}


#' Compare diagonal and correlated-trait OU covariance models
#'
#' Fits a diagonal innovation-covariance null model and a full covariance model
#' to the same shift configuration. A parametric-bootstrap likelihood-ratio test
#' is used because small phylogenetic samples and covariance optimization can
#' make a simple chi-squared reference distribution inaccurate. Zero
#' off-diagonal covariances are interior points of the positive-definite
#' covariance parameter space, not boundary parameters.
#'
#'@param tree ultrametric phylogeny in postorder.
#'@param Y numeric trait matrix aligned to \code{tree$tip.label}.
#'@param shift.configuration fixed shift-edge indices used in both models.
#'@param nboot number of null simulations. Use zero to return fitted models and
#' the observed likelihood ratio without a bootstrap p-value.
#'@param seed optional random seed.
#'@param root.model ancestral-root model.
#'@param alpha.structure alpha structure for the full model. The default allows
#' one alpha per trait, making the diagonal model nested in the full model.
#'@param likelihood.engine likelihood engine for the full model.
#'@param optimizer.starts number of deterministic optimization starts.
#'@param measurement_error logical; estimate additional observation variance.
#'@param input_error optional known observation-error variances.
#'@param quietly suppress progress messages.
#'@param covariance.regularization covariance regularization for the alternative.
#'@param regularization.lambda optional regularization intensity.
#'@param nCores number of bootstrap refits evaluated concurrently where fork
#' parallelism is supported.
#'@param selection whether the covariance bootstrap conditions on the supplied
#' shift configuration or repeats shift discovery in every null simulation.
#'@param search.max.nShifts maximum shifts in a repeated search. By default it
#' is at least one and otherwise the size of \code{shift.configuration}.
#'@param search.strategy search strategy passed to
#' \code{\link{estimate_shift_configuration}} when
#' \code{selection = "full"}.
#'@return An object of class \code{"l1ou_covariance_comparison"} containing the
#' two fits, observed statistic, bootstrap statistics and p-value. Regularized
#' alternatives are explicitly marked as sensitivity analyses.
#'@export
compare_trait_covariance <- function(
        tree, Y, shift.configuration=integer(), nboot=100L, seed=NULL,
        root.model=c("OUfixedRoot", "OUrandomRoot"),
        alpha.structure=c("diagonal", "shared"),
        likelihood.engine=c("auto", "dense", "pruning"),
        optimizer.starts=5L, measurement_error=FALSE, input_error=NULL,
        quietly=TRUE,
        covariance.regularization=c("none", "shrinkage"),
        regularization.lambda=NA_real_, nCores=1L,
        selection=c("conditional", "full"), search.max.nShifts=NULL,
        search.strategy=c("auto", "lasso", "ensemble", "exhaustive")){
    root.model <- match.arg(root.model)
    alpha.structure <- match.arg(alpha.structure)
    likelihood.engine <- match.arg(likelihood.engine)
    covariance.regularization <- match.arg(covariance.regularization)
    selection <- match.arg(selection)
    search.strategy <- match.arg(search.strategy)
    if(is.null(search.max.nShifts)){
        search.max.nShifts <- max(1L, length(shift.configuration))
    }
    search.max.nShifts <- l1ou_integer_argument(
        search.max.nShifts, "search.max.nShifts", 0L
    )
    nCores <- l1ou_integer_argument(nCores, "nCores", 1L)
    quietly <- l1ou_logical_argument(quietly, "quietly")
    validate_l1ou_tree(tree, require.positive.edges=TRUE)
    Y <- as_l1ou_trait_matrix(Y)
    if(!identical(rownames(Y), tree$tip.label)){
        stop("rownames of Y must be identical to tree$tip.label.")
    }
    if(ncol(Y) < 2L) stop("at least two traits are required.")
    null <- fit_OU(
        tree, Y, shift.configuration, criterion="BIC",
        root.model=root.model, measurement_error=measurement_error,
        input_error=input_error, trait.covariance="diagonal"
    )
    alternative <- fit_OU(
        tree, Y, shift.configuration, criterion="BIC",
        root.model=root.model, measurement_error=measurement_error,
        input_error=input_error, trait.covariance="full",
        alpha.structure=alpha.structure, likelihood.engine=likelihood.engine,
        optimizer.starts=optimizer.starts,
        covariance.regularization=covariance.regularization,
        regularization.lambda=regularization.lambda
    )
    statistic <- max(0, 2 * (as.numeric(logLik(alternative)) -
                              as.numeric(logLik(null))))
    nboot <- l1ou_integer_argument(nboot, "nboot", 0L)
    bootstrap.statistics <- numeric(0)
    failures <- 0L
    failure.messages <- character()
    if(nboot > 0L){
        simulations <- simulate(null, nsim=nboot, seed=seed, engine="tree")
        one.replicate <- function(i){
            simulated <- simulations[[i]]
            error.message <- ""
            value <- tryCatch({
                replicate.configuration <- shift.configuration
                if(selection == "full"){
                    selected <- estimate_shift_configuration(
                        tree, simulated,
                        max.nShifts=search.max.nShifts,
                        criterion="BIC",
                        root.model=root.model,
                        search.strategy=search.strategy,
                        measurement_error=measurement_error,
                        input_error=input_error,
                        trait.covariance="diagonal",
                        quietly=TRUE
                    )
                    replicate.configuration <- selected$shift.configuration
                }
                fit.null <- fit_OU(
                    tree, simulated, replicate.configuration, criterion="BIC",
                    root.model=root.model, measurement_error=measurement_error,
                    input_error=input_error, trait.covariance="diagonal"
                )
                fit.alternative <- fit_OU(
                    tree, simulated, replicate.configuration, criterion="BIC",
                    root.model=root.model, measurement_error=measurement_error,
                    input_error=input_error, trait.covariance="full",
                    alpha.structure=alpha.structure,
                    likelihood.engine=likelihood.engine,
                    optimizer.starts=optimizer.starts,
                    compute.hessian=FALSE,
                    covariance.regularization=covariance.regularization,
                    regularization.lambda=regularization.lambda
                )
                max(0, 2 * (as.numeric(logLik(fit.alternative)) -
                             as.numeric(logLik(fit.null))))
            }, error=function(e){
                error.message <<- conditionMessage(e)
                NA_real_
            })
            if(!isTRUE(quietly) && (i %% max(1L, floor(nboot / 10L)) == 0L)){
                message("covariance bootstrap ", i, "/", nboot)
            }
            list(value=value, error=error.message)
        }
        records <- if(nCores > 1L && l1ou_supports_multicore()){
            with_l1ou_thread_limit(1L, l1ou_mclapply(
                seq_along(simulations), one.replicate, mc.cores=nCores
            ))
        } else lapply(seq_along(simulations), one.replicate)
        bootstrap.statistics <- vapply(records, function(x) x$value, numeric(1))
        failure.messages <- vapply(records, function(x) x$error, character(1))
        failures <- sum(!is.finite(bootstrap.statistics))
        bootstrap.statistics <- bootstrap.statistics[is.finite(bootstrap.statistics)]
        if(length(bootstrap.statistics) == 0L){
            stop("all covariance-comparison bootstrap refits failed.")
        }
    }
    exceedances <- sum(bootstrap.statistics >= statistic)
    p.value <- if(nboot > 0L)
        (1 + exceedances) / (1 + length(bootstrap.statistics)) else NA_real_
    p.interval <- if(length(bootstrap.statistics) > 0L){
        stats::binom.test(exceedances, length(bootstrap.statistics))$conf.int
    } else c(NA_real_, NA_real_)
    result <- list(
        null=null,
        alternative=alternative,
        statistic=statistic,
        df=ncol(Y) * (ncol(Y) - 1L) / 2L,
        bootstrap.statistics=bootstrap.statistics,
        p.value=p.value,
        p.value.mcse=if(is.finite(p.value))
            sqrt(p.value * (1 - p.value) /
                 (length(bootstrap.statistics) + 1L)) else NA_real_,
        p.value.interval=p.interval,
        successful=length(bootstrap.statistics),
        failed=failures,
        attempted=nboot,
        selection=selection,
        search.max.nShifts=search.max.nShifts,
        calibrated=identical(covariance.regularization, "none"),
        calibration.note=if(identical(covariance.regularization, "none")) NULL else
            paste0(
                "The alternative was fitted with a covariance penalty; the ",
                "result is a regularized sensitivity comparison, not a ",
                "conventional likelihood-ratio test."
            ),
        failure.messages=sort(
            table(failure.messages[nzchar(failure.messages)]), decreasing=TRUE
        )
    )
    class(result) <- "l1ou_covariance_comparison"
    result
}

#'@export
print.l1ou_covariance_comparison <- function(x, ...){
    cat("Trait-covariance likelihood-ratio comparison\n")
    cat("LR statistic:", signif(x$statistic, 6), "\n")
    if(is.finite(x$p.value)){
        cat("parametric-bootstrap p-value:", signif(x$p.value, 6),
            "(", x$successful, "successful refits)\n")
    }
    invisible(x)
}

#' Summarize uncertainty in bootstrap shift locations
#'
#'@param bootstrap.result result returned by
#' \code{\link{l1ou_bootstrap_support}}.
#'@param tree optional phylogeny used to label edge rows and columns.
#'@return Exact-configuration and tip-partition frequencies, edge inclusion,
#' edge co-selection, and tip co-assignment probabilities.
#'@export
summarize_shift_uncertainty <- function(bootstrap.result, tree=NULL){
    if(!is.list(bootstrap.result)) stop("bootstrap.result must be a list.")
    if(!is.null(tree)) validate_l1ou_tree(tree, require.positive.edges=FALSE)
    shifts <- bootstrap.result$all.shifts
    if(is.null(shifts) || length(shifts) == 0L){
        stop("bootstrap.result contains no successful shift configurations.")
    }
    n.edge <- if(!is.null(tree)) nrow(tree$edge) else{
        max(c(length(bootstrap.result$detection.rate), unlist(shifts), 0L))
    }
    incidence <- matrix(0, nrow=length(shifts), ncol=n.edge)
    for(i in seq_along(shifts)){
        selected <- as.integer(shifts[[i]])
        if(anyNA(selected) || any(selected < 1L | selected > n.edge)){
            stop("bootstrap.result contains an invalid shift-edge index.")
        }
        selected <- unique(selected)
        incidence[i, selected] <- 1
    }
    keys <- vapply(shifts, function(x){
        x <- sort(unique(as.integer(x)))
        if(length(x) == 0L) "none" else paste(x, collapse=",")
    }, character(1))
    counts <- sort(table(keys), decreasing=TRUE)
    configurations <- data.frame(
        configuration=names(counts),
        count=as.integer(counts),
        probability=as.numeric(counts) / length(shifts),
        row.names=NULL, stringsAsFactors=FALSE
    )
    edge.labels <- if(is.null(tree)) as.character(seq_len(n.edge)) else{
        paste0(seq_len(n.edge), ":", tree$edge[, 1L], "->", tree$edge[, 2L])
    }
    colnames(incidence) <- edge.labels
    co.selection <- crossprod(incidence) / nrow(incidence)
    partitions <- if(is.null(tree)) NULL else lapply(
        shifts, function(configuration){
            shift_tip_partition(tree, configuration)
        }
    )
    partition.keys <- if(is.null(tree)) NULL else vapply(
        shifts, function(configuration){
            shift_partition_key(tree, configuration)
        }, character(1)
    )
    partition.frequencies <- NULL
    tip.coassignment <- NULL
    modal.partition <- NULL
    partition.ari <- NULL
    if(!is.null(tree)){
        partition.counts <- sort(table(partition.keys), decreasing=TRUE)
        representative <- match(names(partition.counts), partition.keys)
        partition.frequencies <- data.frame(
            partition=seq_along(partition.counts),
            configuration=vapply(
                shifts[representative],
                function(configuration){
                    if(length(configuration)){
                        paste(sort(configuration), collapse=",")
                    } else "none"
                },
                character(1)
            ),
            count=as.integer(partition.counts),
            probability=as.numeric(partition.counts) / length(shifts),
            key=names(partition.counts),
            row.names=NULL,
            stringsAsFactors=FALSE
        )
        tip.coassignment <- Reduce(`+`, lapply(
            partitions, partition_coassignment_matrix
        )) / length(partitions)
        dimnames(tip.coassignment) <- list(tree$tip.label, tree$tip.label)
        modal.partition <- partitions[[representative[[1L]]]]
        partition.ari <- vapply(
            partitions, adjusted_rand_index, numeric(1), right=modal.partition
        )
    }
    clade.inclusion <- NULL
    if(!is.null(tree)){
        clade.keys <- edge_descendant_tip_keys(tree)$keys
        grouped <- split(seq_len(n.edge), clade.keys)
        clade.inclusion <- vapply(grouped, function(index){
            mean(rowSums(incidence[, index, drop=FALSE]) > 0)
        }, numeric(1))
    }
    result <- list(
        configurations=configurations,
        partitions=partition.frequencies,
        edge.inclusion=stats::setNames(colMeans(incidence), edge.labels),
        co.selection=co.selection,
        tip.coassignment=tip.coassignment,
        modal.partition=modal.partition,
        partition.ari.to.modal=partition.ari,
        clade.inclusion=clade.inclusion,
        number.of.shifts=table(rowSums(incidence)),
        attempted=if(is.null(bootstrap.result$attempted)) length(shifts) else
            bootstrap.result$attempted,
        successful=length(shifts),
        failed=if(is.null(bootstrap.result$failed)) 0L else
            bootstrap.result$failed
    )
    class(result) <- "l1ou_shift_uncertainty"
    result
}

#'@export
print.l1ou_shift_uncertainty <- function(x, ...){
    cat("Bootstrap shift uncertainty\n")
    cat("successful:", x$successful, "/", x$attempted,
        "(failed:", x$failed, ")\n")
    cat("most frequent configurations:\n")
    print(utils::head(x$configurations, 5L), row.names=FALSE)
    if(!is.null(x$partitions)){
        cat("most frequent tip partitions:\n")
        print(utils::head(
            x$partitions[, c("partition", "configuration", "count", "probability")],
            5L
        ), row.names=FALSE)
    }
    invisible(x)
}

#' Information-criterion model averaging over evaluated shift configurations
#'
#'@param model fitted \code{"l1ou"} object containing a search profile.
#'@param delta.max largest score difference retained from the best model.
#'@param max.models optional cap on the number of retained models.
#'@return A list containing normalized weights, refitted models, averaged tip
#' means and optima, averaged process parameters, and edge-inclusion weights.
#'@export
model_average_l1ou <- function(model, delta.max=10, max.models=Inf){
    check_l1ou_object(model)
    if(isTRUE(model$convergent)){
        stop(paste0(
            "model averaging over an unconstrained search profile cannot ",
            "preserve convergent-regime constraints."
        ))
    }
    if(length(delta.max) != 1L || !is.finite(delta.max) || delta.max < 0){
        stop("delta.max must be one finite non-negative number.")
    }
    if(length(max.models) != 1L || !is.numeric(max.models) ||
       is.na(max.models) || (!is.finite(max.models) && max.models != Inf) ||
       (is.finite(max.models) &&
        (max.models < 1 || max.models != floor(max.models)))){
        stop("max.models must be at least one or Inf.")
    }
    profile <- model$profile
    if(is.null(profile) || is.null(profile$scores) ||
       is.null(profile$configurations)){
        stop("the model does not contain a candidate-configuration profile.")
    }
    if(!any(is.finite(profile$scores))){
        stop("the model profile contains no finite score.")
    }
    order <- order(profile$scores)
    scores <- as.numeric(profile$scores[order])
    configurations <- profile$configurations[order]
    keys <- vapply(configurations, function(x){
        paste(sort(as.integer(x)), collapse=",")
    }, character(1))
    unique.keep <- !duplicated(keys)
    scores <- scores[unique.keep]
    configurations <- configurations[unique.keep]
    keep <- is.finite(scores) & scores - min(scores, na.rm=TRUE) <= delta.max
    scores <- scores[keep]
    configurations <- configurations[keep]
    if(is.finite(max.models) && length(scores) > max.models){
        index <- seq_len(as.integer(max.models))
        scores <- scores[index]
        configurations <- configurations[index]
    }
    opt <- make_nested_l1ou_options(model$l1ou.options)
    opt$use.saved.scores <- FALSE
    failure.messages <- character(length(configurations))
    fits <- lapply(seq_along(configurations), function(i){
        configuration <- configurations[[i]]
        tryCatch(
            fit_OU(
                model$tree, model$Y, configuration,
                l1ou.options=opt
            ),
            error=function(e){
                failure.messages[[i]] <<- conditionMessage(e)
                NULL
            }
        )
    })
    successful <- !vapply(fits, is.null, logical(1))
    if(!any(successful)) stop("all candidate-model refits failed.")
    fits <- fits[successful]
    configurations <- configurations[successful]
    profile.scores <- scores[successful]
    scores <- vapply(fits, function(fit) as.numeric(fit$score), numeric(1))
    if(!all(is.finite(scores))){
        stop("one or more successful candidate refits returned non-finite scores.")
    }
    weights <- exp(-0.5 * (scores - min(scores)))
    weights <- weights / sum(weights)
    weighted.array <- function(component){
        Reduce(`+`, Map(function(fit, weight){
            as.matrix(fit[[component]]) * weight
        }, fits, weights))
    }
    n.edge <- nrow(model$tree$edge)
    edge.inclusion <- numeric(n.edge)
    for(i in seq_along(fits)){
        edge.inclusion[as.integer(fits[[i]]$shift.configuration)] <-
            edge.inclusion[as.integer(fits[[i]]$shift.configuration)] + weights[[i]]
    }
    alpha <- Reduce(`+`, Map(function(fit, weight){
        as.numeric(fit$alpha) * weight
    }, fits, weights))
    covariance <- Reduce(`+`, Map(function(fit, weight){
        evolutionary_vcov(fit) * weight
    }, fits, weights))
    partitions <- lapply(configurations, function(configuration){
        shift_tip_partition(model$tree, configuration)
    })
    partition.keys <- vapply(configurations, function(configuration){
        shift_partition_key(model$tree, configuration)
    }, character(1))
    partition.weights <- tapply(weights, partition.keys, sum)
    partition.weights <- sort(partition.weights, decreasing=TRUE)
    tip.coassignment <- Reduce(`+`, Map(function(partition, weight){
        partition_coassignment_matrix(partition) * weight
    }, partitions, weights))
    dimnames(tip.coassignment) <- list(
        model$tree$tip.label, model$tree$tip.label
    )
    result <- list(
        scores=scores,
        profile.scores=profile.scores,
        weights=weights,
        configurations=configurations,
        models=fits,
        mu=weighted.array("mu"),
        optima=weighted.array("optima"),
        alpha=alpha,
        trait.covariance=covariance,
        edge.inclusion=edge.inclusion,
        partition.weights=partition.weights,
        tip.coassignment=tip.coassignment,
        search.diagnostics=model$search.diagnostics,
        failed=sum(!successful),
        failure.messages=sort(
            table(failure.messages[nzchar(failure.messages)]), decreasing=TRUE
        ),
        weight.note=if(isTRUE(model$search.diagnostics$globally.optimal)){
            paste0(
                "Relative exp(-Delta/2) weights over the complete enumerated ",
                "configuration space."
            )
        } else paste0(
            "Relative exp(-Delta/2) weights over the evaluated candidate set; ",
            "they do not account for configurations omitted by the search."
        )
    )
    class(result) <- "l1ou_model_average"
    result
}

#'@export
print.l1ou_model_average <- function(x, ...){
    cat("kfl1ou model average\n")
    cat("models:", length(x$models), "(failed refits:", x$failed, ")\n")
    table <- data.frame(
        score=x$scores,
        weight=x$weights,
        shifts=vapply(x$configurations, length, integer(1))
    )
    print(table, row.names=FALSE)
    invisible(x)
}


#' Repeat shift inference across a phylogenetic tree ensemble
#'
#' Fits the complete shift-search pipeline independently to trees sharing the
#' same extant taxa and aggregates uncertainty by tip partitions, which remain
#' comparable even when edge indices and topologies differ.
#'
#'@param trees non-empty list of phylogenetic trees with the same tip labels.
#'@param Y trait vector or matrix with row names.
#'@param tree.weights optional non-negative weights for the trees.
#'@param ... arguments passed to
#' \code{\link{estimate_shift_configuration}}.
#'@return A tree-ensemble object containing fits, partition weights, mean tip
#' co-assignment probabilities, and pairwise adjusted Rand indices.
#'@export
fit_l1ou_tree_ensemble <- function(trees, Y, tree.weights=NULL, ...){
    if(inherits(trees, "phylo")) trees <- list(trees)
    if(!is.list(trees) || !length(trees)){
        stop("trees must be a non-empty list of phylogenetic trees.")
    }
    Y <- as_l1ou_trait_matrix(Y)
    if(is.null(rownames(Y))) stop("Y must have row names for a tree ensemble.")
    reference <- sort(trees[[1L]]$tip.label)
    if(any(vapply(trees, function(tree){
        !inherits(tree, "phylo") || !identical(sort(tree$tip.label), reference)
    }, logical(1)))){
        stop("all trees must contain exactly the same tip labels.")
    }
    if(!identical(sort(rownames(Y)), reference)){
        stop("rownames(Y) must match the common tree tip labels.")
    }
    if(is.null(tree.weights)) tree.weights <- rep(1, length(trees))
    if(length(tree.weights) != length(trees) || !is.numeric(tree.weights) ||
       any(!is.finite(tree.weights)) || any(tree.weights < 0) ||
       sum(tree.weights) <= 0){
        stop("tree.weights must be finite, non-negative, and have positive sum.")
    }
    tree.weights <- tree.weights / sum(tree.weights)
    errors <- character(length(trees))
    fits <- lapply(seq_along(trees), function(i){
        tree <- reorder(trees[[i]], "postorder")
        response <- Y[tree$tip.label, , drop=FALSE]
        tryCatch(
            estimate_shift_configuration(tree, response, ...),
            error=function(e){
                errors[[i]] <<- conditionMessage(e)
                NULL
            }
        )
    })
    successful <- !vapply(fits, is.null, logical(1))
    if(!any(successful)) stop("shift inference failed on every tree.")
    weights <- tree.weights[successful]
    weights <- weights / sum(weights)
    fits <- fits[successful]
    partitions <- lapply(fits, function(fit){
        partition <- shift_tip_partition(fit$tree, fit$shift.configuration)
        partition[reference]
    })
    keys <- vapply(fits, function(fit){
        shift_partition_key(fit$tree, fit$shift.configuration)
    }, character(1))
    partition.weights <- tapply(weights, keys, sum)
    partition.weights <- sort(partition.weights, decreasing=TRUE)
    tip.coassignment <- Reduce(`+`, Map(function(partition, weight){
        partition_coassignment_matrix(partition) * weight
    }, partitions, weights))
    dimnames(tip.coassignment) <- list(reference, reference)
    pairwise.ari <- outer(seq_along(partitions), seq_along(partitions),
                          Vectorize(function(i, j){
        adjusted_rand_index(partitions[[i]], partitions[[j]])
    }))
    result <- list(
        fits=fits,
        weights=weights,
        partition.weights=partition.weights,
        tip.coassignment=tip.coassignment,
        pairwise.ari=pairwise.ari,
        shift.counts=vapply(fits, function(fit){
            length(fit$shift.configuration)
        }, integer(1)),
        successful=sum(successful),
        failed=sum(!successful),
        failure.messages=sort(table(errors[nzchar(errors)]), decreasing=TRUE)
    )
    class(result) <- "l1ou_tree_ensemble"
    result
}


#'@export
print.l1ou_tree_ensemble <- function(x, ...){
    cat("kfl1ou tree-ensemble sensitivity analysis\n")
    cat("successful:", x$successful, "failed:", x$failed, "\n")
    cat("weighted shift count:", sum(x$weights * x$shift.counts), "\n")
    cat("top tip partitions:\n")
    print(utils::head(x$partition.weights, 5L))
    invisible(x)
}


#' Root-model, alpha-bound, and criterion sensitivity analysis
#'
#' Refits a model over a compact grid of scientifically important assumptions.
#' It can either preserve the selected shifts or repeat the complete search.
#'
#'@param model fitted \code{"l1ou"} object.
#'@param root.models root-state models to compare.
#'@param alpha.upper positive upper bounds to compare.
#'@param criteria information criteria to compare.
#'@param selection condition on the fitted shifts or repeat shift discovery.
#'@param keep.fits retain successful models in an attribute.
#'@return Data frame with scores, likelihoods, shift counts, alpha estimates,
#' and errors for every sensitivity-grid cell.
#'@export
sensitivity_l1ou <- function(
        model,
        root.models=c("OUfixedRoot", "OUrandomRoot"),
        alpha.upper=NULL,
        criteria=NULL,
        selection=c("conditional", "full"),
        keep.fits=FALSE){
    check_l1ou_object(model)
    selection <- match.arg(selection)
    root.models <- match.arg(
        root.models, c("OUfixedRoot", "OUrandomRoot"), several.ok=TRUE
    )
    if(is.null(alpha.upper)){
        fitted.upper <- model$l1ou.options$alpha.upper.bound
        fitted.upper <- fitted.upper[is.finite(fitted.upper) & fitted.upper > 0]
        baseline <- if(length(fitted.upper)) max(fitted.upper) else
            alpha_upper_bound(model$tree)
        alpha.upper <- unique(c(baseline / 2, baseline, baseline * 2))
    }
    if(!is.numeric(alpha.upper) || !length(alpha.upper) ||
       any(!is.finite(alpha.upper)) || any(alpha.upper <= 0)){
        stop("alpha.upper must contain finite positive values.")
    }
    if(is.null(criteria)) criteria <- unique(c(model$l1ou.options$criterion, "BIC"))
    allowed <- if(!is.null(model$trait.covariance)) c("BIC", "pBIC") else
        c("pBIC", "pBICess", "mBIC", "BIC", "AICc")
    if(any(!criteria %in% allowed)){
        stop("one or more criteria are unsupported by this covariance model.")
    }
    grid <- expand.grid(
        root.model=root.models,
        alpha.upper=alpha.upper,
        criterion=criteria,
        stringsAsFactors=FALSE
    )
    fits <- vector("list", nrow(grid))
    errors <- character(nrow(grid))
    for(i in seq_len(nrow(grid))){
        common <- list(
            tree=model$tree,
            Y=model$Y,
            criterion=grid$criterion[[i]],
            root.model=grid$root.model[[i]],
            alpha.upper=grid$alpha.upper[[i]],
            measurement_error=isTRUE(model$l1ou.options$measurement_error),
            input_error=model$l1ou.options$input_error,
            trait.covariance=model$l1ou.options$trait.covariance,
            alpha.structure=model$l1ou.options$alpha.structure,
            covariance.regularization=model$l1ou.options$covariance.regularization,
            regularization.lambda=model$l1ou.options$regularization.lambda,
            likelihood.engine=model$l1ou.options$likelihood.engine,
            optimizer.starts=model$l1ou.options$optimizer.starts,
            compute.hessian=FALSE
        )
        fits[[i]] <- tryCatch({
            if(selection == "full"){
                do.call(
                    estimate_shift_configuration,
                    c(common, list(
                        max.nShifts=model$l1ou.options$max.nShifts,
                        search.strategy=model$l1ou.options$search.strategy,
                        quietly=TRUE
                    ))
                )
            } else{
                do.call(
                    fit_OU,
                    c(common, list(
                        shift.configuration=as.integer(model$shift.configuration)
                    ))
                )
            }
        }, error=function(e){
            errors[[i]] <<- conditionMessage(e)
            NULL
        })
    }
    result <- transform(
        grid,
        score=vapply(fits, function(fit){
            if(is.null(fit)) NA_real_ else as.numeric(fit$score)
        }, numeric(1)),
        logLik=vapply(fits, function(fit){
            if(is.null(fit)) NA_real_ else as.numeric(logLik(fit))
        }, numeric(1)),
        n.shifts=vapply(fits, function(fit){
            if(is.null(fit)) NA_integer_ else length(fit$shift.configuration)
        }, integer(1)),
        alpha=vapply(fits, function(fit){
            if(is.null(fit)) NA_character_ else
                paste(signif(as.numeric(fit$alpha), 6), collapse=",")
        }, character(1)),
        error=errors
    )
    attr(result, "selection") <- selection
    if(isTRUE(keep.fits)) attr(result, "fits") <- fits
    class(result) <- c("l1ou_sensitivity", class(result))
    result
}
