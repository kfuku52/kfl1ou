## generate the W matrix, feature vectors, as described in the doc.
generate_prediction_vec  <-  function(tr, 
                                      shift.configuration, 
                                      conv.regimes, 
                                      alpha, 
                                      ageMatrix=NULL, 
                                      designMatrix=F){

    ## In fact ageMatrix is the approximate design matrix.
    if(is.null(ageMatrix)){
        if ( is.na(alpha) ){
            X   <-  generate_design_matrix(tr, "apprX")
        }else{
            X   <-  generate_design_matrix(tr, "orgX",  alpha = alpha)
        }

        if(designMatrix){
            Cinvh <- t( sqrt_OU_covariance(tr, alpha=alpha, root.model = "OUfixedRoot")$sqrtInvSigma )
            X     <- Cinvh%*%X
        }
        preds <- cbind(1, X[,shift.configuration])
        colnames(preds) <- c(0, shift.configuration)
        Z           <- generate_design_matrix(tr, "simpX")
        template.Z  <- cbind(1, Z[,shift.configuration])
        colnames(template.Z) <- c(0, shift.configuration) 
    }else{
        stopifnot(ncol(ageMatrix)==length(shift.configuration))
        stopifnot(alpha>0)

        preds <- cbind(1, 1-exp(-alpha*ageMatrix))
        colnames(preds) <- c(0, shift.configuration)

        template.Z  <- cbind(1, ageMatrix)
        colnames(template.Z) <- c(0, shift.configuration) 
    }

    ## now the coefficients in the linear regression represent the optimum values.
    ## rather than the shift values.

    for( i in 1:ncol(preds) ){
        for( j in 1:ncol(preds) ){
            if ( i == j )  
                next
            set1 <- which( template.Z[,i] > 0)
            set2 <- which( template.Z[,j] > 0)
            ## if edge j is an ancestor of i
            if ( all(set1 %in% set2) ) 
                preds[ ,j] <- preds[ ,j] - preds[ ,i]
        }
    } 

    W <- numeric()
    if ( length( conv.regimes ) > 0 ){
        for( i in 1:length(conv.regimes) ){
            set1 <- paste(conv.regimes[[i]])
            stopifnot( length(set1) > 0 )

            ## The equality constrain of two optimum values translates 
            ## to combining the corresponding predictors.
            if ( length(set1) > 1){
                W <- cbind( W, rowSums( preds[, paste(set1)] ) ) 
            } else{
                W <- cbind( W,  preds[, paste(set1)] ) 
            }
            colnames(W)[length(W[1,])] <- i
        }
    }

    return(W)
}

normalize_convergent_regimes <- function(regimes, shift.configuration){

    if(is.null(regimes)){
        if(is.null(names(shift.configuration))){
            stop("the convergent regimes must be indicated through the names of the shift.configuration vector.")
        }
        regimes <- list()
        cr.names <- names(shift.configuration)
        idx <- 1L
        for(cr in cr.names){
            regimes[[idx]] <- sort(shift.configuration[which(cr.names == cr)])
            idx <- idx + 1L
        }
    }

    regimes <- lapply(regimes, function(reg){
        sort(unique(as.integer(reg)))
    })
    regimes <- regimes[vapply(regimes, length, integer(1)) > 0L]

    bg.idx <- which(vapply(regimes, function(reg) 0L %in% reg, logical(1)))
    if(length(bg.idx) == 0L){
        regimes <- c(list(0L), regimes)
    } else if(length(bg.idx) > 1L){
        merged.bg <- sort(unique(unlist(regimes[bg.idx])))
        regimes <- c(list(merged.bg), regimes[-bg.idx])
    } else if(bg.idx != 1L){
        regimes <- c(regimes[bg.idx], regimes[-bg.idx])
    }

    return(regimes)
}

map_convergent_regimes_to_trait <- function(regimes, original.shift.configuration,
                                            augmented.shift.configuration){

    original.shift.configuration <- as.integer(original.shift.configuration)
    augmented.shift.configuration <- as.integer(augmented.shift.configuration)

    mapped.regimes <- list()
    seen.keys <- character(0)

    for(reg in regimes){
        reg <- as.integer(reg)
        mapped <- integer(0)

        if(any(reg == 0L)){
            mapped <- c(mapped, 0L)
        }

        edge.reg <- setdiff(reg, 0L)
        if(length(edge.reg) > 0L){
            idx <- match(edge.reg, original.shift.configuration)
            if(any(is.na(idx))){
                stop("convergent regimes do not match with the shift positions.")
            }
            mapped.edges <- augmented.shift.configuration[idx]
            mapped.edges <- mapped.edges[!is.na(mapped.edges)]
            if(length(mapped.edges) > 0L){
                mapped <- c(mapped, mapped.edges)
            }
        }

        mapped <- sort(unique(as.integer(mapped)))
        if(length(mapped) == 0L){
            next
        }

        key <- paste(mapped, collapse = " ")
        if(key %in% seen.keys){
            next
        }

        mapped.regimes[[length(mapped.regimes) + 1L]] <- mapped
        seen.keys <- c(seen.keys, key)
    }

    if(length(mapped.regimes) == 0L){
        return(list(0L))
    }

    bg.idx <- which(vapply(mapped.regimes, function(reg) 0L %in% reg, logical(1)))
    if(length(bg.idx) == 0L){
        mapped.regimes <- c(list(0L), mapped.regimes)
    } else if(length(bg.idx) > 1L){
        merged.bg <- sort(unique(unlist(mapped.regimes[bg.idx])))
        mapped.regimes <- c(list(merged.bg), mapped.regimes[-bg.idx])
    } else if(bg.idx != 1L){
        mapped.regimes <- c(mapped.regimes[bg.idx], mapped.regimes[-bg.idx])
    }

    return(mapped.regimes)
}

dynamic_known_input_error_gls_fit <- function(tree, Y, build_preds,
                                              lower.bound=NA_real_,
                                              upper.bound=NA_real_,
                                              starting.value=NA_real_,
                                              quietly=TRUE,
                                              input_error){

    Y <- as.matrix(Y)
    y <- as.numeric(Y[, 1])
    n <- length(y)
    phy <- reorder(tree, "pruningwise")
    mean.tip.height <- mean(pruningwise.distFromRoot(phy)[seq_len(n)])
    tol <- 1e-10

    if(is.null(names(input_error))){
        names(input_error) <- tree$tip.label
    }
    input_error <- as.numeric(input_error[tree$tip.label])

    fit_for_alpha <- function(alpha){

        preds <- as.matrix(build_preds(alpha))
        d <- ncol(preds)
        if(n <= d){
            return(list(n2llh = .Machine$double.xmax))
        }

        re <- sqrt_OU_covariance(
            tree,
            alpha = alpha,
            root.model = "OUfixedRoot",
            sigma2 = 1,
            input_error = input_error,
            check.order = FALSE,
            check.ultrametric = FALSE
        )
        Linv <- t(re$sqrtInvSigma)
        Xt <- Linv %*% preds
        yt <- drop(Linv %*% y)

        XX <- crossprod(Xt)
        Xy <- crossprod(Xt, yt)

        inv.solve <- tryCatch({
            chol.XX <- chol(XX)
            list(
                invXX = chol2inv(chol.XX),
                betahat = drop(backsolve(chol.XX, forwardsolve(t(chol.XX), Xy)))
            )
        }, error = function(e) {
            invXX <- tryCatch(solve(XX), error = function(e2) NULL)
            if(is.null(invXX)){
                return(NULL)
            }
            list(invXX = invXX, betahat = drop(invXX %*% Xy))
        })

        if(is.null(inv.solve)){
            return(list(n2llh = .Machine$double.xmax))
        }

        fitted.values <- drop(preds %*% inv.solve$betahat)
        residuals <- y - fitted.values
        sigma2hat <- sum((yt - drop(Xt %*% inv.solve$betahat))^2) / n
        if(!is.finite(sigma2hat) || sigma2hat <= 0){
            return(list(n2llh = .Machine$double.xmax))
        }

        logdet <- 2 * as.numeric(determinant(re$sqrtSigma, logarithm = TRUE)$modulus)
        n2llh <- as.numeric(n * log(2 * pi) + n + n * log(sigma2hat) + logdet)
        vcov <- sigma2hat * inv.solve$invXX * n/(n - d)

        list(
            n2llh = n2llh,
            betahat = inv.solve$betahat,
            sigma2hat = sigma2hat,
            vcov = vcov,
            fitted.values = fitted.values,
            residuals = residuals,
            preds = preds,
            d = d
        )
    }

    lower.bound <- ifelse(is.null(lower.bound), NA_real_, as.numeric(lower.bound[[1]]))
    upper.bound <- ifelse(is.null(upper.bound), NA_real_, as.numeric(upper.bound[[1]]))
    starting.value <- ifelse(is.null(starting.value), NA_real_, as.numeric(starting.value[[1]]))

    if(is.na(lower.bound) && is.na(starting.value)){
        lower.opt <- 1e-07 / mean.tip.height
    } else{
        lower.opt <- ifelse(is.na(lower.bound), 0, lower.bound)
    }

    if(is.na(upper.bound) || upper.bound <= 0){
        stop("convergent input_error fits require a strictly positive alpha.upper bound.")
    }

    if(lower.opt == upper.bound){
        alpha.hat <- lower.opt
    } else if(lower.opt > 0){
        opt.res <- optimize(
            function(logalpha) fit_for_alpha(exp(logalpha))$n2llh,
            interval = c(log(lower.opt), log(upper.bound)),
            tol = 1e-6
        )
        alpha.hat <- as.numeric(exp(opt.res$minimum))
    } else{
        opt.res <- optimize(
            function(alpha) fit_for_alpha(alpha)$n2llh,
            interval = c(0, upper.bound),
            tol = 1e-6
        )
        alpha.hat <- as.numeric(opt.res$minimum)
    }

    if(!quietly &&
       (isTRUE(all.equal(alpha.hat, lower.opt, tolerance = tol)) ||
        isTRUE(all.equal(alpha.hat, upper.bound, tolerance = tol)))){
        warning(paste("the estimation of alpha matches the upper/lower bound for this parameter.\n                          You may change the bounds using options \"upper.bound\" and \"lower.bound\".\n"))
    }

    fit <- fit_for_alpha(alpha.hat)
    if(!is.finite(fit$n2llh) || fit$n2llh >= .Machine$double.xmax){
        stop("failed to fit the convergent regime model with the supplied input_error.")
    }

    coefficient.names <- colnames(fit$preds)
    coefficients <- fit$betahat
    names(coefficients) <- coefficient.names
    vcov <- fit$vcov
    colnames(vcov) <- coefficient.names
    rownames(vcov) <- coefficient.names
    p <- fit$d + 2L

    list(
        coefficients = coefficients,
        sigma2 = fit$sigma2hat,
        optpar = alpha.hat,
        sigma2_error = 0,
        logLik = -fit$n2llh/2,
        p = p,
        aic = 2 * p + fit$n2llh,
        vcov = vcov,
        fitted.values = fit$fitted.values,
        residuals = fit$residuals,
        mean.tip.height = mean.tip.height,
        y = y,
        X = fit$preds,
        n = n,
        d = fit$d,
        model = "OUfixedRoot"
    )
}

dynamic_joint_input_measurement_error_gls_fit <- function(tree, Y, build_preds,
                                                          lower.bound=NA_real_,
                                                          upper.bound=NA_real_,
                                                          starting.value=NA_real_,
                                                          quietly=TRUE,
                                                          input_error){

    Y <- as.matrix(Y)
    y <- as.numeric(Y[, 1])
    n <- length(y)
    phy <- reorder(tree, "pruningwise")
    mean.tip.height <- mean(pruningwise.distFromRoot(phy)[seq_len(n)])
    tol <- 1e-10
    objective.ceiling <- .Machine$double.xmax / 1000

    if(is.null(names(input_error))){
        names(input_error) <- tree$tip.label
    }
    input_error <- as.numeric(input_error[tree$tip.label])

    fit_profiled_gls <- function(Sigma, preds){

        Sigma <- 0.5 * (Sigma + t(Sigma))
        chol.Sigma <- tryCatch(chol(Sigma), error = function(e) NULL)
        if(is.null(chol.Sigma)){
            return(NULL)
        }

        Xt <- tryCatch(forwardsolve(t(chol.Sigma), preds), error = function(e) NULL)
        yt <- tryCatch(drop(forwardsolve(t(chol.Sigma), matrix(y, ncol = 1))),
                       error = function(e) NULL)
        if(is.null(Xt) || is.null(yt)){
            return(NULL)
        }

        d <- ncol(preds)
        XX <- crossprod(Xt)
        Xy <- crossprod(Xt, yt)

        inv.solve <- tryCatch({
            chol.XX <- chol(XX)
            list(
                invXX = chol2inv(chol.XX),
                betahat = drop(backsolve(chol.XX, forwardsolve(t(chol.XX), Xy)))
            )
        }, error = function(e) {
            invXX <- tryCatch(solve(XX), error = function(e2) NULL)
            if(is.null(invXX)){
                return(NULL)
            }
            list(invXX = invXX, betahat = drop(invXX %*% Xy))
        })

        if(is.null(inv.solve)){
            return(NULL)
        }

        residuals.whitened <- yt - drop(Xt %*% inv.solve$betahat)
        rss <- sum(residuals.whitened^2)
        if(!is.finite(rss) || rss < 0){
            return(NULL)
        }

        fitted.values <- drop(preds %*% inv.solve$betahat)
        residuals <- y - fitted.values
        logdet <- 2 * sum(log(diag(chol.Sigma)))
        n2llh <- as.numeric(n * log(2 * pi) + logdet + rss)
        vcov.scale <- if(n > d) rss/(n - d) else 1

        list(
            n2llh = n2llh,
            betahat = inv.solve$betahat,
            vcov = inv.solve$invXX * vcov.scale,
            fitted.values = fitted.values,
            residuals = residuals,
            preds = preds,
            d = d
        )
    }

    lower.bound <- ifelse(is.null(lower.bound), NA_real_, as.numeric(lower.bound[[1]]))
    upper.bound <- ifelse(is.null(upper.bound), NA_real_, as.numeric(upper.bound[[1]]))
    starting.value <- ifelse(is.null(starting.value), NA_real_, as.numeric(starting.value[[1]]))

    if(is.na(lower.bound) && is.na(starting.value)){
        alpha.lower <- 1e-07 / mean.tip.height
    } else{
        alpha.lower <- ifelse(is.na(lower.bound), 0, lower.bound)
    }
    alpha.upper <- upper.bound
    if(is.na(alpha.upper) || alpha.upper <= 0){
        stop("convergent measurement_error fits require a strictly positive alpha.upper bound.")
    }

    alpha.hat <- NA_real_
    optimize.alpha <- !isTRUE(all.equal(alpha.lower, alpha.upper, tolerance = tol))
    alpha.log.scale <- optimize.alpha && alpha.lower > 0
    if(optimize.alpha){
        alpha.start <- ifelse(is.na(starting.value), max(0.5/mean.tip.height, alpha.lower), starting.value)
        alpha.start <- min(max(alpha.start, alpha.lower), alpha.upper)
    } else{
        alpha.hat <- alpha.lower
    }

    covariance.cache <- new.env(parent = emptyenv())

    get_phylo_covariance <- function(alpha){

        key <- paste0("alpha:", format(alpha, digits = 16, scientific = FALSE))
        if(exists(key, envir = covariance.cache, inherits = FALSE)){
            return(get(key, envir = covariance.cache, inherits = FALSE))
        }

        re <- sqrt_OU_covariance(
            tree,
            alpha = alpha,
            root.model = "OUfixedRoot",
            sigma2 = 1,
            check.order = FALSE,
            check.ultrametric = FALSE
        )
        Sigma <- tcrossprod(re$sqrtSigma)
        Sigma <- 0.5 * (Sigma + t(Sigma))
        assign(key, Sigma, envir = covariance.cache)
        Sigma
    }

    baseline <- try(
        dynamic_known_input_error_gls_fit(
            tree,
            Y,
            build_preds = build_preds,
            lower.bound = lower.bound,
            upper.bound = upper.bound,
            starting.value = starting.value,
            quietly = quietly,
            input_error = input_error
        ),
        silent = TRUE
    )

    sigma2.start <- if(inherits(baseline, "try-error")){
        max(stats::var(y), 1e-08)
    } else{
        max(as.numeric(baseline$sigma2), 1e-08)
    }
    sigma2_error.start <- max(median(input_error, na.rm = TRUE) * 0.25, 0)

    decode_parameters <- function(par){

        idx <- 0L
        if(optimize.alpha){
            idx <- idx + 1L
            alpha <- if(alpha.log.scale) exp(par[[idx]]) else par[[idx]]
        } else{
            alpha <- alpha.hat
        }

        idx <- idx + 1L
        sigma2 <- exp(par[[idx]])
        idx <- idx + 1L
        sigma2_error <- par[[idx]]

        list(alpha = alpha, sigma2 = sigma2, sigma2_error = sigma2_error)
    }

    fit_for_parameters <- function(alpha, sigma2, sigma2_error){

        if(!is.finite(sigma2) || sigma2 <= 0 || !is.finite(sigma2_error) || sigma2_error < 0){
            return(list(n2llh = objective.ceiling))
        }

        preds <- as.matrix(build_preds(alpha))
        d <- ncol(preds)
        if(n <= d){
            return(list(n2llh = objective.ceiling))
        }

        phylo.cov <- get_phylo_covariance(alpha)
        total.cov <- sigma2 * phylo.cov
        diag(total.cov) <- diag(total.cov) + input_error + sigma2_error

        fit <- fit_profiled_gls(total.cov, preds)
        if(is.null(fit)){
            return(list(n2llh = objective.ceiling))
        }

        fit$sigma2 <- sigma2
        fit$sigma2_error <- sigma2_error
        fit
    }

    objective <- function(par){
        prm <- decode_parameters(par)
        fit_for_parameters(prm$alpha, prm$sigma2, prm$sigma2_error)$n2llh
    }

    par <- lower <- upper <- numeric()
    if(optimize.alpha){
        if(alpha.log.scale){
            par <- c(par, log(alpha.start))
            lower <- c(lower, log(alpha.lower))
            upper <- c(upper, log(alpha.upper))
        } else{
            par <- c(par, alpha.start)
            lower <- c(lower, alpha.lower)
            upper <- c(upper, alpha.upper)
        }
    }

    par <- c(par, log(sigma2.start), sigma2_error.start)
    lower <- c(lower, log(.Machine$double.eps), 0)
    upper <- c(upper, Inf, Inf)

    opt.res <- optim(
        par,
        fn = objective,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper
    )

    prm <- decode_parameters(opt.res$par)
    if(is.na(alpha.hat)){
        alpha.hat <- prm$alpha
    }
    fit <- fit_for_parameters(prm$alpha, prm$sigma2, prm$sigma2_error)

    if(!is.finite(fit$n2llh) || fit$n2llh >= objective.ceiling){
        stop("failed to fit the convergent regime model with the supplied input_error and measurement_error.")
    }

    if(!quietly && optimize.alpha &&
       (isTRUE(all.equal(alpha.hat, alpha.lower, tolerance = tol)) ||
        isTRUE(all.equal(alpha.hat, alpha.upper, tolerance = tol)))){
        warning(paste("the estimation of alpha matches the upper/lower bound for this parameter.\n                          You may change the bounds using options \"upper.bound\" and \"lower.bound\".\n"))
    }

    coefficient.names <- colnames(fit$preds)
    coefficients <- fit$betahat
    names(coefficients) <- coefficient.names
    vcov <- fit$vcov
    colnames(vcov) <- coefficient.names
    rownames(vcov) <- coefficient.names
    p <- fit$d + 3L

    list(
        coefficients = coefficients,
        sigma2 = fit$sigma2,
        optpar = alpha.hat,
        sigma2_error = fit$sigma2_error,
        logLik = -fit$n2llh/2,
        p = p,
        aic = 2 * p + fit$n2llh,
        vcov = vcov,
        fitted.values = fit$fitted.values,
        residuals = fit$residuals,
        mean.tip.height = mean.tip.height,
        y = y,
        X = fit$preds,
        n = n,
        d = fit$d,
        model = "OUfixedRoot"
    )
}

convergent_error_interface_CR <- function(tr, Y, conv.regimes = list(), alpha=NA,
                                          fixed.alpha=FALSE, opt,
                                          shift.configuration = opt$shift.configuration,
                                          input_error = opt$input_error){

    fixed.alpha <- isTRUE(fixed.alpha) || isTRUE(opt$fixed.alpha)

    build_preds <- function(alpha.value){
        preds <- generate_prediction_vec(
            tr,
            shift.configuration,
            conv.regimes,
            alpha = alpha.value,
            ageMatrix = NULL
        )
        if(fixed.alpha){
            preds <- ifelse(preds > 0, 1, 0)
            colnames(preds) <- colnames(generate_prediction_vec(
                tr,
                shift.configuration,
                conv.regimes,
                alpha = alpha.value,
                ageMatrix = NULL
            ))
        }
        preds
    }

    lower.bound <- if(fixed.alpha) alpha else alpha/100
    upper.bound <- if(fixed.alpha) alpha else max(alpha, opt$alpha.upper.bound + .Machine$double.eps)
    starting.value <- alpha

    if(isTRUE(opt$measurement_error)){
        if(is.null(input_error)){
            input_error <- stats::setNames(rep(0, length(tr$tip.label)), tr$tip.label)
        }
        return(dynamic_joint_input_measurement_error_gls_fit(
            tr,
            Y,
            build_preds = build_preds,
            lower.bound = lower.bound,
            upper.bound = upper.bound,
            starting.value = starting.value,
            quietly = TRUE,
            input_error = input_error
        ))
    }

    return(dynamic_known_input_error_gls_fit(
        tr,
        Y,
        build_preds = build_preds,
        lower.bound = lower.bound,
        upper.bound = upper.bound,
        starting.value = starting.value,
        quietly = TRUE,
        input_error = input_error
    ))
}

phylolm_interface_CR  <-  function(tr, Y, conv.regimes = list(), alpha=NA, fixed.alpha=FALSE, opt,
                                   shift.configuration = opt$shift.configuration,
                                   input_error = opt$input_error){

    if( isTRUE(opt$measurement_error) || !is.null(input_error) ){
        return(convergent_error_interface_CR(
            tr,
            Y,
            conv.regimes = conv.regimes,
            alpha = alpha,
            fixed.alpha = fixed.alpha,
            opt = opt,
            shift.configuration = shift.configuration,
            input_error = input_error
        ))
    }

    preds <- generate_prediction_vec(tr, shift.configuration, conv.regimes, alpha, ageMatrix=NULL)
    prev.val <- options()$warn
    options(warn = -1) # dangerous!

    if( is.null(opt$fixed.alpha) ) 
	    opt$fixed.alpha  <- FALSE;


    if(fixed.alpha || opt$fixed.alpha){
	    preds <- ifelse(preds>0,1,0)
	    fit <-  phylolm(Y~preds-1,
			    phy  = tr,
			    model = "OUfixedRoot",
			    starting.value = list(alpha=alpha),
			    lower.bound = alpha, 
			    upper.bound = alpha)

    }else{
	    fit <-  phylolm_CR(Y~preds-1, 
			       phy  = tr, 
			       model = "OUfixedRoot",
			       sc=shift.configuration,
			       cr=conv.regimes,
			       starting.value=alpha,
			       upper.bound=max(alpha, opt$alpha.upper.bound + .Machine$double.eps),
			       lower.bound=alpha/100
			       )
    }
    options(warn = prev.val)
    return(fit)
}

## compute the AICc score
cmp_AICc_CR  <-  function(tree, Y, conv.regimes, alpha, opt){

    shift.configuration <- opt$shift.configuration
    conv.regimes <- normalize_convergent_regimes(conv.regimes, shift.configuration)
    stopifnot( length(alpha) == ncol(Y) )

    nShifts    <- length( shift.configuration )
    nShiftVals <- length( conv.regimes ) -1## conv.regimes has intercept as an optimum value
    nTips      <- length( tree$tip.label )

    p <- nShifts + (nShiftVals + 3 + extra_error_df(opt))*ncol(Y)
    N <- nTips*ncol(Y)
    df.1 <- 2*p + (2*p*(p+1))/(N-p-1) 
    if( p > N-2)  ##  for this criterion we should have p < N.
        return(Inf)
    df.2 <- 0
    score <- df.1
    for( i in 1:ncol(Y)){
        r <- get_data(tree, Y, shift.configuration, opt, i)
        trait.regimes <- map_convergent_regimes_to_trait(
            conv.regimes,
            shift.configuration,
            r$augmented.s.c
        )
        fit <- phylolm_interface_CR(
            r$tr,
            r$y,
            trait.regimes,
            alpha = alpha[[i]],
            opt = opt,
            shift.configuration = r$s.c,
            input_error = r$input_error
        )
        if ( all( is.na( fit) ) ){ return(Inf) } 
        score <- score  -2*fit$logLik + df.2
    }
    return(score)
}

## compute the BIC score
cmp_BIC_CR <- function(tree, Y, conv.regimes, alpha, opt){

    shift.configuration <- opt$shift.configuration
    conv.regimes <- normalize_convergent_regimes(conv.regimes, shift.configuration)
    stopifnot( length(alpha) == ncol(Y) )

    nEdges     <- Nedge(tree)
    nTips      <- length(tree$tip.label)
    nShifts    <- length(shift.configuration)
    nShiftVals <- length( conv.regimes ) - 1 
    nVariables <- ncol(Y)

    df.1  <- log(nTips)*(nShiftVals)
    score <- df.1
    #alpha <- sigma2 <- logLik <- rep(0, nVariables)

    for( i in 1:nVariables ){

        df.2 <- log(nTips)*(nShifts + 3 + extra_error_df(opt))
        r <- get_data(tree, Y, shift.configuration, opt, i)
        trait.regimes <- map_convergent_regimes_to_trait(
            conv.regimes,
            shift.configuration,
            r$augmented.s.c
        )
        fit  <- phylolm_interface_CR(
            r$tr,
            r$y,
            trait.regimes,
            alpha = alpha[[i]],
            opt = opt,
            shift.configuration = r$s.c,
            input_error = r$input_error
        )
        if ( all(is.na(fit)) ){ return(Inf) } 
        score <- score  -2*fit$logLik + df.2
    }
    return( score )
}



## compute the pBIC score
cmp_pBIC_CR  <-  function(tree, Y, conv.regimes, alpha, opt){

    shift.configuration <- opt$shift.configuration
    conv.regimes <- normalize_convergent_regimes(conv.regimes, shift.configuration)
    nShifts = length(shift.configuration)
    nEdges  = Nedge(tree)
    nTips   = length(tree$tip.label)

    df.1   <- 0
    df.1   <- 2*(nShifts)*log(nEdges-1)
    score  <- df.1
    #alpha  <- sigma2 <- logLik <- rep(0, ncol(Y))

    for(i in 1:ncol(Y)){
        r <- get_data(tree, Y, shift.configuration, opt, i)
        trait.regimes <- map_convergent_regimes_to_trait(
            conv.regimes,
            shift.configuration,
            r$augmented.s.c
        )
        fit <- phylolm_interface_CR(
            r$tr,
            r$y,
            trait.regimes,
            alpha = alpha[[i]],
            opt = opt,
            shift.configuration = r$s.c,
            input_error = r$input_error
        )
        fit2 <- phylolm_interface_CR(
            r$tr,
            r$y,
            trait.regimes,
            alpha = alpha[[i]],
            fixed.alpha = TRUE,
            opt = opt,
            shift.configuration = r$s.c,
            input_error = r$input_error
        )
        if( all( is.na(fit) ) ){
           return(Inf)
        } 
        varY  <- var(as.numeric(r$y))
        ld    <- as.numeric(determinant(fit2$vcov * (fit$n - fit$d)/(varY*fit$n),
                                        logarithm = TRUE)$modulus)
        df.2  <- (2 + extra_error_df(opt))*log(nrow(r$y)) - ld
        score <- score  -2*fit$logLik + df.2
    }
    return( score )
}


cmp_model_score_CR <- function(tree, Y, regimes=NULL, alpha=NA, opt){

    shift.configuration <- opt$shift.configuration
    regimes <- normalize_convergent_regimes(regimes, shift.configuration)

    if( opt$criterion == "AICc"){
        score <- cmp_AICc_CR(tree, Y, conv.regimes = regimes, alpha=alpha, opt=opt)
    } else if( opt$criterion == "pBIC"){
        score <- cmp_pBIC_CR(tree, Y, conv.regimes = regimes, alpha=alpha, opt=opt)
    } else if( opt$criterion == "BIC"){
        score <- cmp_BIC_CR(tree, Y, conv.regimes = regimes, alpha=alpha,  opt=opt)
    } else
        stop("undefined criterion for convergent evolution!")

    return(score)
}


generate_relation  <- function(tr, shift.configuration){

    s.p     = shift.configuration
    n.s.p   = length(s.p)
    nEdges  = Nedge(tr)

    M  <- numeric()
    tmp.s.p <- s.p
    for ( s1 in s.p ){
        tmp.s.p <- setdiff(tmp.s.p, s1)
        #for ( s2 in setdiff(s.p, s1) )
        for ( s2 in tmp.s.p )
        {
            nr <- rep(0, nEdges) 
            ## B M ADDed
            nr[[s1]] <-  1
            nr[[s2]] <- -1
            nr       <- nr[s.p]

            ## removing rows with only one +-1, we assume that previous step took care of redundancy.
            if ( length( which( abs(nr)>0) ) < 2)
                next

            M <- rbind(M, nr)
            rownames(M)[[length(rownames(M))]] <- paste0(s1," ",s2)
        }
    }

    colnames(M) <- s.p
    stopifnot( nrow(M) > 0 )
    return( M )
}

find_convergent_regimes  <-  function(tr, Y, alpha, criterion, regimes){
    stopifnot(ncol(Y)==1)
    stopifnot(all( row.names(Y) == tr$tip.label))
    l1ou_require_genlasso()

    #alpha <- eModel$alpha
    Cinvh   <- t( sqrt_OU_covariance(tr, alpha=alpha, root.model = "OUfixedRoot")$sqrtInvSigma )
    #Cinvh   <- t( cmp.OU.covariance(tr, alpha=alpha)$D ) 
    Y  <- Cinvh%*%Y

    shift.configuration <- unlist(regimes)
    X   <-  generate_prediction_vec(tr, shift.configuration, alpha=alpha, conv.regimes=regimes, designMatrix=TRUE)
    #X   <-  X[,-1]
    X   <- cbind(X,1)

    M   <- generate_relation(tr, 1:length(regimes))
    M   <- cbind(M,0)

    ###I multipled YY by a number to scale it up. If I don't genlasso doesn't return the whole solution path :S
    out <- NULL
    last.error <- NULL
    genlasso_fit <- getFromNamespace("genlasso", "genlasso")
    for( eps in c(0, 1e-8, 1e-6, 1e-4) ){
        out <- tryCatch(
            genlasso_fit(100*Y, X, D=M, svd=T, eps=eps, approx=F, verbose=F),
            error = function(e) {
                last.error <<- e
                NULL
            }
        )
        if( !is.null(out) ){
            break
        }
        msg <- conditionMessage(last.error)
        if( !grepl("eps must be positive|column rank deficient", msg) ){
            stop(last.error)
        }
    }
    if( is.null(out) ){
        stop(last.error)
    }

    ### adding lambda=0 and removing the intercept
    out$beta   <- cbind(out$beta, coef(out, lambda=0)$beta )
    out$lambda <- c(out$lambda, 0)

    rownames(out$beta) <- colnames(X)

    ## removing the intercpet value from coefficients and adding it to a new vector.
    spots         <- length( out$beta[,1] )
    out$intercept <- out$beta[    spots, ]
    out$beta      <- out$beta[   -spots, , drop=FALSE]
    M             <- M       [ , -spots, drop=FALSE ]

    out$M <- M
    return(out)
}

## This method finds the convergent evolution model that maximizes the 
## criterion such as AICc. To find the optimum, it combines shifts into a 
## convergent regime or splits a regime if that decreases the criterion until 
## no progress.  NOTE: CR is a short for convergent regime.
estimate_convergent_regimes_surface  <-  function(model, opt){

    criterion <- opt$criterion
    model$l1ou.options$criterion = criterion # model is returned. for correct print and summary
    Y         <- as.matrix(model$Y)
    tr        <- model$tree

    sc.prev <- sc  <- model$shift.configuration
    min.regimes    <- as.list(c(0,sc))
    min.score      <- cmp_model_score_CR(tr, Y, regimes=min.regimes, alpha=model$alpha, opt=opt)
    prev.min.score <- min.score
    ## elist represents the edgelist format of the regimes graph.
    ## At the beginning each regime forms a vertex with a self-loop. 
    elist.ref  <-  numeric()
    for(u in c(0,sc)){ elist.ref <- rbind( elist.ref, as.character(c(u,u)) ) }
    current.num.cc <- nrow(elist.ref)

    for( iter in 1:(2*length(sc)) ){

        has.progress <- FALSE
        elist.min <- elist.ref
        ##NOTE: merge, add the edge (u,v), regimes if the IC decreases the most.

        run.list <- list()
	list.idx <- 1

        for( u in c(0, sc) ){ 
            for( v in c(0, sc) ){

                ##NOTE: test if we can add (u, v)
                if( u == v ){ next }
                ##FIXME: The convergence of the immediate shifts after the root
                ##to the background is pointless so don't check them.

                ## add the edge (u,v) to the graph. 
                elist <- as.matrix( rbind(elist.ref, as.character(c(u,v)) ) )
                cc    <- connected_components_from_edgelist(elist)
                ## check if it connects two connected components
                ## of the graph. If not, then it is redundant.
                if( length(cc) >= current.num.cc ){ next }
                ## extract the connected components 
                regimes <-  lapply(cc, as.numeric)

                ## name each cr as the smallest shift index in it 
                sc.tmp  <-  sort(sc)
                for( thelist in lapply(regimes, sort) ){
                    names(sc.tmp)[sc.tmp %in% thelist]  <-  thelist[[1]] 
                }
                if(identical(names(sc.tmp), names(sc.prev))){ next }
                sc.prev <- sc.tmp

		if(!opt$parallel.computing){
			score   <-  cmp_model_score_CR(tr, Y, regimes, model$alpha, opt=opt)
			if( min.score > score ){
				min.score    <- score
				min.regimes  <- regimes
				elist.min    <- elist
				has.progress <- TRUE
			}
		}else{
			run.list[[list.idx]] <- list(elist=elist, regimes=regimes)
			list.idx <- list.idx+1
		}


            }
        }

	if( length(run.list)>0 && opt$parallel.computing ){

		RE.list <- l1ou_mclapply( run.list, FUN=function(X){
					    return ( cmp_model_score_CR(tr, Y, X$regimes, model$alpha, opt=opt) )
			       }, mc.cores=opt$nCores)

		for(idx in 1:length(RE.list) ){
			IN <- run.list[[idx]]
			score   <- RE.list[[idx]] 

			if( min.score > score ){
				min.score    <- score
				min.regimes  <- IN$regimes
				elist.min    <- IN$elist
				has.progress <- TRUE
			}
		}
	}

        if( !has.progress ){
            break
        }

        current.num.cc <- length(min.regimes) 
        elist.ref      <- elist.min

        ## break a CR if it increases the score
        ##FIXME: we don't need to check the very last edge we added
        if( has.progress){
            for( e.idx in 1:length(elist.ref[,1]) ){

                u <- elist.ref[e.idx, 1]
                v <- elist.ref[e.idx, 2]
                if(u==v){next}

                elist <- elist.ref[ -e.idx, ]
                cc    <- connected_components_from_edgelist(elist)

                if( length(cc) <= current.num.cc ){ next }

                regimes <- lapply(cc, as.numeric)
                score   <- cmp_model_score_CR(tr, Y, regimes, model$alpha, opt=opt)

                if( min.score > score ){
                    min.score   <- score
                    elist.min   <- elist
                    min.regimes <- regimes
                }
            }
        }

        elist.ref      <- elist.min
        current.num.cc <- length(min.regimes) 

        #if no progress then terminate
        if( !has.progress ){ break }
        prev.min.score <- min.score
    }
    
    sc <- model$shift.configuration
    counter <- 1
    for (reg in min.regimes) {
        ## 0 represents the background(intercept) 
        background <- (0 %in% reg)
        for (item in sort(reg)) {
            names(sc)[which(sc == item)] <- ifelse(background, 0, counter)
        }
        counter <- counter + ifelse(background, 0, 1)
    }

    model$shift.configuration <- sc
    model$score               <- min.score 
    return(model)
}



#' Detects convergent regimes under an OU model
#'
#' Takes a model previously estimated by \code{\link{estimate_shift_configuration}},
#' including one or more traits and a configuration of evolutionary shifts, and detect which of these regime shifts
#' are convergent.
#'
#'@param model fitted object of class l1ou returned by \pkg{kfl1ou}.
#'@param criterion information criterion for model selection (see Details in \code{\link{configuration_ic}}).
#'@param method search method for finding convergent regimes. ``rr'' is based on
#'  the optionally installed \code{genlasso} package for regularized linear regression.
#'  Currenly, this method can only accept a single trait.
#'  The default ``backward'' method is a heuristic similar to \code{surface_backward}
#'  in the \code{surface} package,
#'  using backward steps to repeatedly merge similar regimes into convergent regimes.
#'  Models fitted with \code{measurement_error} or \code{input_error} are
#'  currently supported only by \code{method = "backward"}.
#'@param fixed.alpha indicates if the alpha parameters should be optimized during likelihood fitting.
#'@param nCores maximum total CPU budget for \code{kfl1ou}. If \code{nCores=1}
#' then it runs sequentially. Otherwise, when fork-based parallelism is
#' available, \code{kfl1ou} may use up to \code{nCores} forked workers via
#' \code{mclapply}, while BLAS/OpenMP threads are limited when supported so the
#' overall computation respects this budget when possible.
#'
#'@details \code{nCores} follows the same total-budget semantics as in
#'         \code{\link{estimate_shift_configuration}}. When process parallelism
#'         is used, worker processes run with BLAS/OpenMP threads limited to 1
#'         when supported, and the previous thread settings are restored on exit.
#'         Models created by \pkg{kfl1ou} retain class \code{"l1ou"} for
#'         backward compatibility with upstream code.
#'
#'@examples
#' 
#'data("lizard.traits", "lizard.tree")
#'keep <- lizard.tree$tip.label[1:15]
#'tree <- drop.tip(lizard.tree, setdiff(lizard.tree$tip.label, keep))
#'tree <- reorder(tree, "postorder")
#'Y <- as.matrix(lizard.traits[keep, 1, drop = FALSE])
#' ## first fit a model to find individual shifts (no convergence assumed):
#'fit_ind <- estimate_shift_configuration(tree, Y, criterion="AICc", max.nShifts=2)
#' ## then detect which of these shifts are convergent:
#'fit_conv <- estimate_convergent_regimes(fit_ind, criterion="AICc")
#'plot(fit_conv)
#'
#'@seealso   \code{\link{estimate_shift_configuration}}
#'
#'@export
estimate_convergent_regimes  <-  function(model, 
                                        criterion=c("AICc", "pBIC", "BIC"),
                                        method=c("backward", "rr"),
					fixed.alpha=FALSE,
					nCores=1
                                     ){
    opt <- list()
    opt$method <- match.arg(method)
    opt$criterion <- match.arg(criterion)
    model$l1ou.options$criterion = criterion
    opt$shift.configuration <- model$shift.configuration
    opt$alpha.upper.bound  <- model$l1ou.options$alpha.upper.bound
    opt$measurement_error <- isTRUE(model$l1ou.options$measurement_error)
    opt$input_error <- model$l1ou.options$input_error
    opt$multivariate.missing <- isTRUE(model$l1ou.options$multivariate.missing)
    opt$tree.list <- model$l1ou.options$tree.list
    opt$quietly <- TRUE
    opt$fixed.alpha <- fixed.alpha

    opt$nCores <- nCores
    opt$parallel.computing <- FALSE
    if( opt$nCores > 1){
        if(!l1ou_supports_multicore()){
            warning("fork-based parallel execution is unavailable; running sequentially.", immediate=TRUE)
            opt$nCores <- 1
        }else{
	    opt$parallel.computing <- TRUE
        }
    }
    opt$ageMatrix <- generate_design_matrix(model$tree, "apprX")[, model$shift.configuration]

    thread.limit <- resolve_l1ou_thread_limit(opt$nCores, opt$parallel.computing)
    return(with_l1ou_thread_limit(thread.limit, {
        if(opt$method == "backward"){
            return(estimate_convergent_regimes_surface(model, opt=opt))
        }

        if(isTRUE(opt$measurement_error) || !is.null(opt$input_error)){
            stop("estimate_convergent_regimes(method=\"rr\") does not yet support measurement_error or input_error.")
        }


        ## doing this for genlasso. FIXME: change the solver as it seems unstable.
        Y   <-  32*model$Y/lnorm(model$Y,l=2)
        Y   <-  as.matrix(Y)
        tr  <-  model$tree

        stopifnot( ncol(Y) == 1 ) # this method only works for univariate trait
        if( length(model$shift.configuration) < 2 ){
            model$score <- cmp_model_score_CR(tr, model$Y,
                                              regimes = as.list(c(0, model$shift.configuration)),
                                              alpha = model$alpha, opt=opt)
            return(model)
        }

        c.regimes <- prev.regimes <- all.regimes <- list()
        c.regimes[1:length(model$shift.configuration)] <- model$shift.configuration
        min.cr.regimes <- c.regimes
        min.score <- cmp_model_score_CR(tr, Y, c(list(0), c.regimes), model$alpha, opt=opt)
        prev.min.score <- min.score
        ar.counter <- 1
        
        ## similar to "backward" method. But here we may combine several shifts into convergent regimes
        ## at the same time therefore it is faster.
        for(iter in 1:length(model$shift.configuration) ){
            out  <-  find_convergent_regimes(tr, Y, model$alpha, opt$criterion, regimes = c.regimes )
            for(num.digits in c(12,13,15,16)){
                for( idx in 1:length(out$beta[1,]) ){
    
                new.est  <- round( out$M%*%out$beta[,idx], digits = num.digits)
                conv.reg <- rownames(out$M)[ which( new.est == 0) ]
    
                elist    <- numeric()
                for(e in conv.reg){## converting to numerical matrix.
                    elist <- rbind( elist, unlist(strsplit(e," ") ) )
                }
                for(e in paste(1:length(c.regimes)) ){
                    elist <- rbind( elist, c(e,e))
                }
    
                ## extracting connected components as the convergent regimes 
                cc      <- connected_components_from_edgelist(elist)
                s.p.tmp <- model$shift.configuration
                regimes <- list()
                counter <- 1

                for(i in 1:length(cc) ){

                    regimes[[counter]] <- numeric()
                    the.cc <- as.numeric(cc[[i]])
                    for( vv in the.cc ){
                        regimes[[counter]] <- c(regimes[[counter]], c.regimes[[vv]])
                    }
                    s.p.tmp <- setdiff( s.p.tmp, regimes[[counter]] )
                    counter <- counter + 1
                }
    
                ## s.p.tmp contains shifts that are not in the graph
                for( s in s.p.tmp){
                    regimes[[counter]] <- s
                    counter <- counter + 1
                }
    
                if( identical(prev.regimes, regimes) )
                    next
                prev.regimes <- regimes
    
                res <- lapply( X=all.regimes, FUN=function(x){ return(identical(regimes,x)) }  )
                if( any( unlist(res)) )
                    next
    
                all.regimes[[ ar.counter ]] <- regimes
                ar.counter <- ar.counter + 1
    
                score <- cmp_model_score_CR(tr, Y, c(list(0), regimes), model$alpha, opt=opt)
    
                if( min.score > score ){
                    min.score      <- score
                    min.cr.regimes <- regimes
                    min.digits     <- num.digits
                }
            }
        }

        c.regimes <- min.cr.regimes
        if( min.score == prev.min.score )
            break;
        prev.min.score <- min.score

        }

        sc <- model$shift.configuration
        counter <- 1
        for( reg in c.regimes ){
            for( item in sort(reg) ){
                names(sc)[which(sc==item)] <- counter 
            }
            counter <- counter + 1
        }

        model$shift.configuration <- sc
        model$score  <-  cmp_model_score_CR(tr, model$Y, c(list(0), c.regimes), model$alpha, opt=opt)

        return(model)
    }))
}






 
