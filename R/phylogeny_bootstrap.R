#'
#' Computes bootstrap support for shift positions
#'
#' Takes a given shift configuration previously detected from data along with shift magnitudes
#' and OU parameters, to calculate bootstrap support for shift positions. 
#' The non-parametric bootstrap procedure calculates phylogenetically-uncorrelated standardized residuals,
#' one at each node. These residuals are sampled with replacement, then mapped back onto the tree
#' to create bootstrap replicates. Each replicate is analyzed with
#' \code{kfl1ou} and user-specified options.
#'
#'@param model an object output by \code{\link{estimate_shift_configuration}}. 
#'@param nItrs number of independent iterations (bootstrap replicates).
#'@param multicore logical. If TRUE, bootstrap replicates are distributed over
#' parallel workers.
#'@param nCores maximum total CPU budget for the bootstrap run. If
#' \code{multicore=TRUE}, at most \code{nCores} bootstrap replicates are
#' processed in parallel; each worker evaluates \code{kfl1ou} sequentially with
#' BLAS/OpenMP threads limited to 1 when supported, so nested \code{kfl1ou} calls
#' do not oversubscribe the machine.
#'@param quietly logical. If FALSE, a summary of each iteration will be printed out.
#'
#'
#'@return A list with \code{detection.rate}, a vector of edge-wise detection
#' proportions, and \code{all.shifts}, the shift configurations retained from
#' successful bootstrap replicates.
#'
#'
#'@details The results of sequential and parallel runs are not necessarily equal,
#'         because different seeds might be used for different bootstrap
#'         replicates. Fork-based parallel execution is used only on supported
#'         platforms. To change options for the analysis of
#'         each bootstrap replicate, like the information criterion or the
#'         maximum allowed number of shifts, modify \code{model$l1ou.options}.
#'         This component keeps its upstream name for backward compatibility.
#'         Replicates that fail during model fitting are skipped; the function
#'         throws an error only if all bootstrap replicates fail.
#'
#'
#'@examples
#' 
#' data(lizard.traits, lizard.tree)
#' keep <- lizard.tree$tip.label[1:15]
#' tree <- drop.tip(lizard.tree, setdiff(lizard.tree$tip.label, keep))
#' tree <- reorder(tree, "postorder")
#' Y <- lizard.traits[keep, 1]
#' eModel <- estimate_shift_configuration(tree, Y, criterion="AICc", max.nShifts=1)
#' result <- l1ou_bootstrap_support(eModel, nItrs=1)
#' # using only 1 replicate is vastly insufficient in general,
#' # but used here to make the illustrative example run faster.
#' nEdges <- Nedge(tree)
#' e.w <- rep(1,nEdges) 
#' if (length(eModel$shift.configuration) > 0) {
#'   e.w[eModel$shift.configuration] <- 3
#' }
#' e.l <- round(result$detection.rate*100, digits=1)
#' # to avoid annotating edges with support at or below 10%
#' e.l <- ifelse(e.l>10, paste0(e.l,"%"), NA)
#' plot(eModel, edge.label=e.l, edge.ann.cex=0.7, edge.label.ann=TRUE, cex=0.5, label.offset=0.02, edge.width=e.w)
#'
#'
#' Y <- lizard.traits[keep, 1:2]
#' eModel <- estimate_shift_configuration(tree, Y, criterion="AICc", max.nShifts=1)
#' result <- l1ou_bootstrap_support(eModel, nItrs=1, multicore=TRUE, nCores=2)
#' result$detection.rate
#'
#'@seealso   \code{\link{estimate_shift_configuration}}
#'
#'@export
l1ou_bootstrap_support <- function(model, nItrs=100, multicore=FALSE, nCores = 2, quietly=TRUE){
 
    if (!inherits(model, "l1ou"))  stop("object \"model\" is not of class \"l1ou\".")

    if(multicore && !l1ou_supports_multicore()){
        warning("fork-based parallel execution is unavailable; running sequentially.", immediate.=TRUE)
        multicore = FALSE
    }

    tree = model$tree
    if(ncol(model$Y)==1){
        return(bootstrap_support_univariate(tree=tree, model=model, nItrs=nItrs, 
                                            multicore=multicore, nCores=nCores, quietly=quietly))
    }
    if(ncol(model$Y)>1){
        return(bootstrap_support_multivariate(tree=tree, model=model, nItrs=nItrs, 
                                              multicore=multicore, nCores=nCores, quietly=quietly))
    }
}

is_failed_bootstrap_fit <- function(x){
    return(length(x) == 1L && isTRUE(is.na(x)))
}

finalize_bootstrap_results <- function(all.shift.configurations, detection.vec){
    valid.count <- 0L
    keep <- rep(TRUE, length(all.shift.configurations))

    for(i in seq_along(all.shift.configurations)){
        shifts <- all.shift.configurations[[i]]
        if(is_failed_bootstrap_fit(shifts)){
            keep[[i]] <- FALSE
            next
        }

        valid.count <- valid.count + 1L
        detection.vec[shifts] <- detection.vec[shifts] + 1
    }

    if(valid.count == 0L){
        stop("all bootstrap replicates failed.")
    }

    return(list(
        detection.rate = detection.vec / valid.count,
        all.shifts = all.shift.configurations[keep]
    ))
}

bootstrap_support_univariate <- function(tree, model, nItrs, multicore=FALSE, nCores=2, quietly=FALSE){

    input.error <- get_trait_input_error(model$l1ou.options, 1, tree=tree)
    RE    = sqrt_OU_covariance(tree, alpha=model$alpha, 
                               root.model = model$l1ou.options$root.model,
                               sigma2 = model$sigma2,
                               sigma2_error = model$sigma2_error,
                               input_error = input.error,
                               check.order=F, check.ultrametric=F)

    C.IH  = t(RE$sqrtInvSigma)
    C.H   = RE$sqrtSigma

    Y     = model$Y
    YY    = C.IH%*%(Y - model$mu )

    seed.vec <- sample(.Machine$integer.max, nItrs+1, replace=TRUE)

    detection.vec = rep(0, nrow(tree$edge))
    all.shift.configurations <- list()
    nested.opt <- make_nested_l1ou_options(model$l1ou.options)

    if(quietly==FALSE)
        print(paste0("iteration #:nShifts:shift configuraitons"))

    valid.count <- 0
    if(multicore == FALSE){
        for(itr in 1:nItrs){
            set.seed(seed.vec[[itr]])
            YYstar = sample(YY, replace = TRUE)
            Ystar  = as.matrix( (C.H%*%YYstar) + model$mu )
            rownames(Ystar) <- rownames(Y)

            eM  <-  tryCatch({
                estimate_shift_configuration(tree, Ystar, l1ou.options =model$l1ou.options)
            }, error = function(e) {
                if(!quietly)
                    message("l1OU error, return NA")
                return(NA) }  )
            if(is_failed_bootstrap_fit(eM)) {next}

            valid.count <- valid.count + 1

            detection.vec[eM$shift.configuration] = detection.vec[eM$shift.configuration] + 1
            all.shift.configurations[[length(all.shift.configurations) + 1L]] <- eM$shift.configuration

            if(quietly==FALSE){
                print(paste0("iteration ", itr, ":", length(eM$shift.configuration),":", 
                             paste0(eM$shift.configuration, collapse=" ") ) )
            }

        }
        if(valid.count == 0L)
            stop("all bootstrap replicates failed.")
        set.seed(seed.vec[[nItrs+1]])
        return(list( detection.rate=(detection.vec/valid.count), all.shifts=all.shift.configurations))
    }


    all.shift.configurations = with_l1ou_thread_limit(1L, 
        l1ou_mclapply(X=1:nItrs, FUN=function(itr){

                     set.seed(seed.vec[[itr]])
                     YYstar = sample(YY, replace = TRUE)
                     Ystar  = as.matrix( (C.H%*%YYstar) + model$mu )
                     rownames(Ystar) <- rownames(Y)

                     eM  <-  tryCatch({
                         estimate_shift_configuration(tree, Ystar, l1ou.options = nested.opt)
                     }, error = function(e) {
                         if(!quietly)
                             message("l1OU error, return NA")
                         return(NA) }  )

                     if(is_failed_bootstrap_fit(eM)) {return(NA)}

                     if(quietly==FALSE){
                         print(paste0("iteration ", itr, ":", length(eM$shift.configuration),":", 
                                      paste0(eM$shift.configuration, collapse=" ") ) )
                     }
                     return(eM$shift.configuration)
           }, mc.cores = nCores))

    set.seed(seed.vec[[nItrs+1]])
    return(finalize_bootstrap_results(all.shift.configurations, detection.vec))
}

bootstrap_support_multivariate <- function(tree, model, nItrs, multicore=FALSE, nCores=2, quietly=FALSE){

    Y = as.matrix(model$Y)
    stopifnot( length(model$alpha) == ncol(Y) )

    seed.vec <- sample(.Machine$integer.max, nItrs+1, replace=TRUE)

    YY        = Y
    C.Hlist   = list()
    for( idx in 1:ncol(Y) ){
        input.error <- get_trait_input_error(model$l1ou.options, idx, tree=tree)
        RE    = sqrt_OU_covariance(tree, alpha = model$alpha[[idx]], 
                                   root.model = model$l1ou.options$root.model,
                                   sigma2 = model$sigma2[[idx]],
                                   sigma2_error = model$sigma2_error[[idx]],
                                   input_error = input.error,
                                   check.order=F, check.ultrametric=F) 
        C.IH  = t(RE$sqrtInvSigma) 
        C.Hlist[[idx]] = RE$sqrtSigma
        YY[, idx]      = C.IH%*%(Y[, idx] - model$mu[ ,idx])

    }
    if(quietly==FALSE)
        print(paste0("iteration #:nShifts:shift configuraitons"))

    detection.vec = rep(0, nrow(tree$edge))
    all.shift.configurations <- list()
    nested.opt <- make_nested_l1ou_options(model$l1ou.options)

    valid.count <- 0
    if( multicore == FALSE ){
        for(itr in 1:nItrs){

            set.seed(seed.vec[[itr]])
            Ystar   = YY
            idx.vec = sample(1:nrow(YY), replace = TRUE)
            for( idx in 1:ncol(YY) ){
                YYstar        = YY[idx.vec, idx]
                Ystar[, idx]  = (C.Hlist[[idx]] %*% YYstar) + model$mu[, idx] 
            }
            rownames(Ystar) <- rownames(Y)

            eM  <-  tryCatch({
                estimate_shift_configuration(tree, Ystar,  l1ou.options=model$l1ou.options)
            }, error = function(e) {
                if(!quietly)
                    message("l1OU error, return NA")
                return(NA) }  )

            if(is_failed_bootstrap_fit(eM)) {next}

            if(quietly==FALSE){
                print(paste0("iteration ", itr, ":", length(eM$shift.configuration),":", 
                             paste0(eM$shift.configuration, collapse=" ") ) )
            }

            valid.count  <- valid.count + 1
            detection.vec[eM$shift.configuration] = detection.vec[eM$shift.configuration] + 1
            all.shift.configurations[[length(all.shift.configurations) + 1L]] <- eM$shift.configuration
        }
        if(valid.count == 0L)
            stop("all bootstrap replicates failed.")
        set.seed(seed.vec[[nItrs+1]])
        return(list( detection.rate=(detection.vec/valid.count), all.shifts=all.shift.configurations))
    }

    all.shift.configurations = with_l1ou_thread_limit(1L,
        l1ou_mclapply(X=1:nItrs, FUN=function(itr){
                     Ystar   = YY
                     set.seed(seed.vec[[itr]])
                     idx.vec = sample(1:nrow(YY), replace = TRUE)
                     for( idx in 1:ncol(YY) ){
                         YYstar        = YY[idx.vec, idx]
                         Ystar[, idx]  = (C.Hlist[[idx]] %*% YYstar) + model$mu[, idx] 
                     }
                     eM  <-  tryCatch({
                         estimate_shift_configuration(tree, Ystar, l1ou.options = nested.opt)
                     }, error = function(e) {
                         if(!quietly)
                             message("l1OU error, return NA")
                         return(NA) }  )

                     if(is_failed_bootstrap_fit(eM)) {return(NA)}

                     if(quietly==FALSE){
                         print(paste0("iteration ", itr, ":", length(eM$shift.configuration),":", 
                                      paste0(eM$shift.configuration, collapse=" ") ) )
                     }

                     return(eM$shift.configuration)
        }, mc.cores = nCores))

    set.seed(seed.vec[[nItrs+1]]) ## To make sure after both mclapply and for-loop we have same seed for the reproducibility  
    return(finalize_bootstrap_results(all.shift.configurations, detection.vec))
}
