#' Compare diagonal and correlated-trait OU covariance models
#'
#' Fits a diagonal innovation-covariance null model and a full covariance model
#' to the same shift configuration. A parametric-bootstrap likelihood-ratio test
#' is used because zero correlations place the null on a constrained part of the
#' covariance parameter space and small phylogenetic samples can invalidate a
#' simple chi-squared reference distribution.
#'
#'@param tree ultrametric phylogeny in postorder.
#'@param Y numeric trait matrix aligned to code{tree$tip.label}.
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
#'@return An object of class code{"l1ou_covariance_comparison"} containing the
#' two fits, observed statistic, bootstrap statistics and p-value.
#'@export
compare_trait_covariance <- function(
        tree, Y, shift.configuration=integer(), nboot=100L, seed=NULL,
        root.model=c("OUfixedRoot", "OUrandomRoot"),
        alpha.structure=c("diagonal", "shared"),
        likelihood.engine=c("auto", "dense", "pruning"),
        optimizer.starts=5L, measurement_error=FALSE, input_error=NULL,
        quietly=TRUE){
    root.model <- match.arg(root.model)
    alpha.structure <- match.arg(alpha.structure)
    likelihood.engine <- match.arg(likelihood.engine)
    Y <- as.matrix(Y)
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
        optimizer.starts=optimizer.starts
    )
    statistic <- max(0, 2 * (as.numeric(logLik(alternative)) -
                              as.numeric(logLik(null))))
    nboot <- as.integer(nboot[[1L]])
    if(is.na(nboot) || nboot < 0L) stop("nboot must be a non-negative integer.")
    bootstrap.statistics <- numeric(0)
    failures <- 0L
    if(nboot > 0L){
        simulations <- simulate(null, nsim=nboot, seed=seed)
        bootstrap.statistics <- vapply(seq_along(simulations), function(i){
            simulated <- simulations[[i]]
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
                    compute.hessian=FALSE
                )
                max(0, 2 * (as.numeric(logLik(fit.alternative)) -
                             as.numeric(logLik(fit.null))))
            }, error=function(e) NA_real_)
            if(!isTRUE(quietly) && (i %% max(1L, floor(nboot / 10L)) == 0L)){
                message("covariance bootstrap ", i, "/", nboot)
            }
            value
        }, numeric(1))
        failures <- sum(!is.finite(bootstrap.statistics))
        bootstrap.statistics <- bootstrap.statistics[is.finite(bootstrap.statistics)]
        if(length(bootstrap.statistics) == 0L){
            stop("all covariance-comparison bootstrap refits failed.")
        }
    }
    result <- list(
        null=null,
        alternative=alternative,
        statistic=statistic,
        df=ncol(Y) * (ncol(Y) - 1L) / 2L,
        bootstrap.statistics=bootstrap.statistics,
        p.value=if(nboot > 0L)
            (1 + sum(bootstrap.statistics >= statistic)) /
                (1 + length(bootstrap.statistics)) else NA_real_,
        successful=length(bootstrap.statistics),
        failed=failures
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
#' code{link{l1ou_bootstrap_support}}.
#'@param tree optional phylogeny used to label edge rows and columns.
#'@return Configuration frequencies, edge inclusion probabilities, and an
#' edge co-selection probability matrix.
#'@export
summarize_shift_uncertainty <- function(bootstrap.result, tree=NULL){
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
        selected <- selected[selected >= 1L & selected <= n.edge]
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
    list(
        configurations=configurations,
        edge.inclusion=stats::setNames(colMeans(incidence), edge.labels),
        co.selection=co.selection,
        number.of.shifts=table(rowSums(incidence)),
        successful=length(shifts)
    )
}

#' Information-criterion model averaging over evaluated shift configurations
#'
#'@param model fitted code{"l1ou"} object containing a search profile.
#'@param delta.max largest score difference retained from the best model.
#'@param max.models optional cap on the number of retained models.
#'@return A list containing normalized weights, refitted models, averaged tip
#' means and optima, averaged process parameters, and edge-inclusion weights.
#'@export
model_average_l1ou <- function(model, delta.max=10, max.models=Inf){
    check_l1ou_object(model)
    profile <- model$profile
    if(is.null(profile) || is.null(profile$scores) ||
       is.null(profile$configurations)){
        stop("the model does not contain a candidate-configuration profile.")
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
    weights <- exp(-0.5 * (scores - min(scores)))
    weights <- weights / sum(weights)
    opt <- make_nested_l1ou_options(model$l1ou.options)
    opt$use.saved.scores <- FALSE
    fits <- lapply(configurations, function(configuration){
        tryCatch(
            fit_OU(
                model$tree, model$Y, configuration,
                l1ou.options=opt
            ),
            error=function(e) NULL
        )
    })
    successful <- !vapply(fits, is.null, logical(1))
    if(!any(successful)) stop("all candidate-model refits failed.")
    fits <- fits[successful]
    configurations <- configurations[successful]
    scores <- scores[successful]
    weights <- weights[successful] / sum(weights[successful])
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
        vcov(fit) * weight
    }, fits, weights))
    result <- list(
        scores=scores,
        weights=weights,
        configurations=configurations,
        models=fits,
        mu=weighted.array("mu"),
        optima=weighted.array("optima"),
        alpha=alpha,
        trait.covariance=covariance,
        edge.inclusion=edge.inclusion
    )
    class(result) <- "l1ou_model_average"
    result
}
