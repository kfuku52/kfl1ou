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
        regularization.lambda=NA_real_, nCores=1L){
    root.model <- match.arg(root.model)
    alpha.structure <- match.arg(alpha.structure)
    likelihood.engine <- match.arg(likelihood.engine)
    covariance.regularization <- match.arg(covariance.regularization)
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
                fit.null <- fit_OU(
                    tree, simulated, shift.configuration, criterion="BIC",
                    root.model=root.model, measurement_error=measurement_error,
                    input_error=input_error, trait.covariance="diagonal"
                )
                fit.alternative <- fit_OU(
                    tree, simulated, shift.configuration, criterion="BIC",
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
#'@return Configuration frequencies, edge inclusion probabilities, and an
#' edge co-selection probability matrix.
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
        edge.inclusion=stats::setNames(colMeans(incidence), edge.labels),
        co.selection=co.selection,
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
        failed=sum(!successful),
        failure.messages=sort(
            table(failure.messages[nzchar(failure.messages)]), decreasing=TRUE
        ),
        weight.note=paste0(
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
