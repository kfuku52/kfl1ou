#' Fit kfl1ou from replicated species measurements
#'
#' Aggregates replicate observations to species means and propagates the
#' estimated sampling variance of each mean as known tip-specific input error.
#' This separates within-species measurement/sampling variation from the
#' evolutionary process instead of estimating both from one value per tip.
#'
#'@param tree phylogenetic tree.
#'@param Y replicate-by-trait numeric matrix.
#'@param species character vector mapping each replicate row to a tip label.
#'@param shift.configuration optional fixed shift configuration. \code{NULL}
#' performs shift discovery.
#'@param singleton.error fallback variance of a mean for species-trait cells
#' with one observation. \code{"pooled"} uses the pooled within-species
#' variance divided by the cell count; a non-negative number is used directly.
#'@param ... arguments passed to
#' \code{\link{estimate_shift_configuration}} or \code{\link{fit_OU}}.
#' Do not supply \code{input_error}.
#'@return A fitted \code{"l1ou"} object with replicate summaries attached.
#'@export
fit_l1ou_replicates <- function(
        tree, Y, species, shift.configuration=NULL,
        singleton.error=c("pooled", "zero"), ...){
    validate_l1ou_tree(tree, require.positive.edges=TRUE)
    Y <- as.matrix(Y)
    storage.mode(Y) <- "double"
    if(nrow(Y) != length(species)){
        stop("species must contain one tip label per row of Y.")
    }
    species <- as.character(species)
    if(anyNA(species) || any(!nzchar(species))){
        stop("species contains missing or empty labels.")
    }
    if(any(!species %in% tree$tip.label)){
        stop("species contains labels absent from tree$tip.label.")
    }
    missing.tips <- setdiff(tree$tip.label, species)
    if(length(missing.tips)){
        stop("every tree tip needs at least one replicate observation.")
    }
    if(is.null(colnames(Y))) colnames(Y) <- paste0("trait", seq_len(ncol(Y)))
    if(is.character(singleton.error)){
        singleton.error <- match.arg(singleton.error)
    } else if(length(singleton.error) != 1L || !is.numeric(singleton.error) ||
              !is.finite(singleton.error) || singleton.error < 0){
        stop('singleton.error must be "pooled", "zero", or a non-negative number.')
    }

    means <- matrix(
        NA_real_, nrow=length(tree$tip.label), ncol=ncol(Y),
        dimnames=list(tree$tip.label, colnames(Y))
    )
    counts <- input.error <- means
    pooled.variance <- vapply(seq_len(ncol(Y)), function(j){
        centered <- unlist(lapply(split(Y[, j], species), function(values){
            values <- values[is.finite(values)]
            if(length(values) < 2L) numeric() else values - mean(values)
        }), use.names=FALSE)
        if(length(centered) < 2L) 0 else sum(centered^2) /
            max(1L, length(centered) - length(unique(species)))
    }, numeric(1))
    for(i in seq_along(tree$tip.label)){
        rows <- species == tree$tip.label[[i]]
        for(j in seq_len(ncol(Y))){
            values <- Y[rows, j]
            values <- values[is.finite(values)]
            counts[i, j] <- length(values)
            if(!length(values)) next
            means[i, j] <- mean(values)
            if(length(values) >= 2L){
                input.error[i, j] <- stats::var(values) / length(values)
            } else if(is.numeric(singleton.error)){
                input.error[i, j] <- singleton.error
            } else if(singleton.error == "zero"){
                input.error[i, j] <- 0
            } else{
                input.error[i, j] <- pooled.variance[[j]]
            }
        }
    }
    if(any(rowSums(is.finite(means)) == 0L)){
        stop("at least one tip has no finite trait observations.")
    }
    input.error[!is.finite(input.error)] <- 0
    call.arguments <- list(...)
    if("input_error" %in% names(call.arguments)){
        stop("input_error is computed from replicates and must not be supplied in ...")
    }
    fit <- if(is.null(shift.configuration)){
        do.call(
            estimate_shift_configuration,
            c(list(tree=tree, Y=means, input_error=input.error), call.arguments)
        )
    } else{
        do.call(
            fit_OU,
            c(list(
                tree=tree, Y=means,
                shift.configuration=shift.configuration,
                input_error=input.error
            ), call.arguments)
        )
    }
    fit$replicate.summary <- list(
        means=means,
        counts=counts,
        input.error=input.error,
        pooled.within.species.variance=stats::setNames(
            pooled.variance, colnames(Y)
        )
    )
    fit
}


#' Measurement-error sensitivity profile
#'
#' Replaces an estimated shared observation variance by a grid of fixed known
#' variances, refits the selected mean model, and reports the likelihood shape.
#'
#'@param model fitted \code{"l1ou"} model.
#'@param multipliers non-negative multipliers of the fitted shared error
#' variance. Include zero to compare against no additional error.
#'@return Data frame of multiplier, log likelihood, score, and refit status.
#'@export
profile_measurement_error_l1ou <- function(
        model, multipliers=c(0, 0.25, 0.5, 1, 2, 4)){
    check_l1ou_object(model)
    if(is.null(model$sigma2_error)){
        stop("model does not contain an estimated measurement-error variance.")
    }
    if(!is.numeric(multipliers) || !length(multipliers) ||
       any(!is.finite(multipliers)) || any(multipliers < 0)){
        stop("multipliers must be finite non-negative numbers.")
    }
    Y <- as.matrix(model$Y)
    base.error <- normalize_input_error(
        model$tree, Y, model$l1ou.options$input_error
    )
    if(is.null(base.error)) base.error <- matrix(
        0, nrow=nrow(Y), ncol=ncol(Y), dimnames=dimnames(Y)
    )
    sigma.error <- rep(as.numeric(model$sigma2_error), length.out=ncol(Y))
    fits <- lapply(multipliers, function(multiplier){
        fixed.error <- sweep(
            base.error, 2L, multiplier * sigma.error, `+`
        )
        opt <- make_nested_l1ou_options(model$l1ou.options)
        opt$measurement_error <- FALSE
        opt$input_error <- fixed.error
        tryCatch(
            fit_OU(
                model$tree, model$Y,
                as.integer(model$shift.configuration),
                l1ou.options=opt
            ),
            error=function(e) e
        )
    })
    data.frame(
        multiplier=multipliers,
        logLik=vapply(fits, function(fit){
            if(inherits(fit, "error")) NA_real_ else as.numeric(logLik(fit))
        }, numeric(1)),
        score=vapply(fits, function(fit){
            if(inherits(fit, "error")) NA_real_ else as.numeric(fit$score)
        }, numeric(1)),
        error=vapply(fits, function(fit){
            if(inherits(fit, "error")) conditionMessage(fit) else ""
        }, character(1)),
        row.names=NULL
    )
}
