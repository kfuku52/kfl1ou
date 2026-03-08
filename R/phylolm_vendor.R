# Vendored and adapted from the CRAN package phylolm 2.6.5 (GPL >= 2).
# Original package authors are credited in NOTICE.
# Modified for kfl1ou on 2026-03-08.

phylolm <- function(formula, data = list(), phy,
                    model = c("BM", "OUrandomRoot", "OUfixedRoot",
                              "lambda", "kappa", "delta", "EB", "trend"),
                    lower.bound = NULL, upper.bound = NULL,
                    starting.value = NULL, measurement_error = FALSE,
                    boot = 0, full.matrix = TRUE, save = FALSE,
                    REML = FALSE, ...) {

    if(boot != 0){
        stop("boot > 0 is not supported by the vendored phylolm implementation.")
    }
    if(isTRUE(full.matrix) || isTRUE(save)){
        invisible(NULL)
    }

    if(!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
    model <- match.arg(model)
    if((model == "trend") && is.ultrametric(phy)){
        stop("the trend is unidentifiable for ultrametric trees.")
    }
    if((model == "lambda") && measurement_error){
        stop("the lambda transformation and measurement error cannot be used together: they are not distinguishable")
    }
    if(is.null(phy$edge.length)) stop("the tree has no branch lengths.")
    if(is.null(phy$tip.label)) stop("the tree has no tip labels.")
    tol <- 1e-10

    mf <- model.frame(formula = formula, data = data)
    if(is.null(rownames(mf))){
        if(nrow(mf) != length(phy$tip.label)){
            stop("number of rows in the data does not match the number of tips in the tree.")
        }
        warning("the data has no names, order assumed to be the same as tip labels in the tree.\n")
    } else {
        taxa_without_data <- setdiff(phy$tip.label, rownames(mf))
        if(length(taxa_without_data) > 0){
            warning("will drop from the tree ", length(taxa_without_data), " taxa with missing data")
            phy <- drop.tip(phy, taxa_without_data)
        }
        if(length(phy$tip.label) < 2){
            stop("names of the data do not match with tip labels.")
        }
        taxa_notin_tree <- setdiff(rownames(mf), phy$tip.label)
        if(length(taxa_notin_tree) > 0){
            warning(length(taxa_notin_tree), " taxa not in the tree: their data will be ignored")
            mf <- mf[-which(rownames(mf) %in% taxa_notin_tree), , drop = FALSE]
        }
        ordr <- match(phy$tip.label, rownames(mf))
        if(any(is.na(ordr))){
            stop("data names do not match with the tip labels.\n")
        }
        mf <- mf[ordr, , drop = FALSE]
    }

    X <- model.matrix(attr(mf, "terms"), data = mf)
    y <- model.response(mf)
    d <- ncol(X)
    if(is.matrix(y) && ncol(y) > 1){
        stop(paste0("The response vector y in the formula is multivariate (it has several columns).\n",
                    "  Please fit each column one by one: 'phylolm' can only handle a simple (univariate) response vector y."))
    }

    phy <- reorder(phy, "pruningwise")
    n <- length(phy$tip.label)
    N <- nrow(phy$edge)
    ROOT <- n + 1L
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    externalEdge <- des <= n

    OU <- c("OUrandomRoot", "OUfixedRoot")
    flag <- 0
    D <- NULL

    if(model %in% OU){
        D <- numeric(n)
        if(!is.ultrametric(phy)){
            flag <- 1
            dis <- pruningwise.distFromRoot(phy)
            Tmax <- max(dis[1:n])
            D <- Tmax - dis[1:n]
            D <- D - mean(D)
            phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] +
                D[des[externalEdge]]
            Tmax <- Tmax + min(D)
        }
    }

    if(model == "trend"){
        trend <- pruningwise.distFromRoot(phy)[1:n]
        X <- cbind(X, trend)
        d <- d + 1
    }

    dis <- pruningwise.distFromRoot(phy)[1:n]
    Tmax <- mean(dis)

    bounds.default <- matrix(c(1e-7 / Tmax, 50 / Tmax,
                               1e-7, 1,
                               1e-6, 1,
                               1e-5, 3,
                               -3 / Tmax, 0,
                               1e-16, 1e16),
                             ncol = 2, byrow = TRUE)
    rownames(bounds.default) <- c("alpha", "lambda", "kappa", "delta", "rate",
                                  "sigma2_error")
    colnames(bounds.default) <- c("min", "max")

    starting.values.default <- c(0.5 / Tmax, 0.5, 0.5, 0.5, -1 / Tmax, 1)
    names(starting.values.default) <- c("alpha", "lambda", "kappa", "delta",
                                        "rate", "sigma2_error")

    get_value_param <- function(values, values.default, param) {
        if(length(values) == 1 && is.null(names(values))){
            if(param == "sigma2_error"){
                val <- NULL
            } else {
                val <- values
            }
        } else {
            val <- values[[param]]
        }
        if(is.null(val)) return(values.default[[param]])
        return(val)
    }

    get_value_model <- function(values, values.default, model, measurement_error) {
        if(model != "BM"){
            if(model %in% OU) param <- "alpha"
            if(model == "lambda") param <- "lambda"
            if(model == "kappa") param <- "kappa"
            if(model == "delta") param <- "delta"
            if(model == "EB") param <- "rate"
            vals <- get_value_param(values, values.default, param)
        } else {
            vals <- NULL
        }
        if(measurement_error){
            vals <- c(vals, get_value_param(values, values.default, "sigma2_error"))
        }
        return(vals)
    }

    lower.bound <- get_value_model(lower.bound, bounds.default[, 1], model, measurement_error)
    upper.bound <- get_value_model(upper.bound, bounds.default[, 2], model, measurement_error)
    starting.value <- get_value_model(starting.value, starting.values.default, model, measurement_error)

    prm <- list(myname = starting.value[1])
    names(prm) <- model
    if(model %in% OU) names(prm) <- "alpha"
    if(model == "EB") names(prm) <- "rate"
    if(measurement_error){
        if(model %in% c("BM", "trend")){
            names(prm) <- "sigma2_error"
        } else {
            prm[["sigma2_error"]] <- starting.value[2]
        }
    }

    ole <- 4 + 2 * d + d * d
    loglik <- function(parameters, y, X) {
        tree <- transf.branch.lengths(phy, model,
                                      parameters = parameters,
                                      check.pruningwise = FALSE,
                                      check.ultrametric = FALSE,
                                      D = D,
                                      check.names = FALSE)$tree
        if(flag){
            y <- exp(-parameters$alpha * D) * y
            X <- exp(-parameters$alpha * D) * X
        }

        tmp <- .C("threepoint_l1ou",
                  as.integer(N), as.integer(n), as.integer(phy$Nnode),
                  as.integer(1), as.integer(d), as.integer(ROOT),
                  as.double(tree$root.edge), as.double(tree$edge.length),
                  as.integer(des), as.integer(anc),
                  as.double(as.vector(y)), as.double(as.vector(X)),
                  result = double(ole),
                  PACKAGE = "kfl1ou")$result

        comp <- list(
            vec11 = tmp[2],
            y1 = tmp[3],
            yy = tmp[4],
            X1 = tmp[5:(4 + d)],
            XX = matrix(tmp[(5 + d):(ole - d)], d, d),
            Xy = tmp[(ole - d + 1):ole],
            logd = tmp[1]
        )
        invXX <- solve(comp$XX)
        betahat <- invXX %*% comp$Xy
        sigma2hat <- as.numeric((comp$yy - 2 * t(betahat) %*% comp$Xy +
            t(betahat) %*% comp$XX %*% betahat) / n)
        if(sigma2hat < 0){
            resdl <- X %*% betahat - y
            tmpyy <- .C("threepoint_l1ou",
                        as.integer(N), as.integer(n), as.integer(phy$Nnode),
                        as.integer(1), as.integer(d), as.integer(ROOT),
                        as.double(tree$root.edge), as.double(tree$edge.length),
                        as.integer(des), as.integer(anc),
                        as.double(as.vector(resdl)), as.double(as.vector(X)),
                        result = double(ole),
                        PACKAGE = "kfl1ou")$result[4]
            sigma2hat <- tmpyy / n
        }
        vcov <- sigma2hat * invXX * n / (n - d)
        if(!REML){
            n2llh <- as.numeric(n * log(2 * pi) + n + n * log(sigma2hat) + comp$logd)
        } else {
            sigma2hat <- sigma2hat * n / (n - d)
            ldXX <- determinant(comp$XX, logarithm = TRUE)$modulus
            n2llh <- as.numeric((n - d) * (log(2 * pi) + 1 + log(sigma2hat)) +
                comp$logd + ldXX)
        }
        if(flag){
            n2llh <- n2llh + parameters$alpha * 2 * sum(D)
        }
        return(list(n2llh = n2llh,
                    betahat = as.vector(betahat),
                    sigma2hat = sigma2hat,
                    vcov = vcov))
    }

    lower <- lower.bound
    upper <- upper.bound
    start <- starting.value

    if((model %in% c("BM", "trend")) && (!measurement_error)){
        BMest <- loglik(prm, y, X)
        results <- list(coefficients = BMest$betahat,
                        sigma2 = BMest$sigma2hat,
                        optpar = NULL,
                        sigma2_error = 0,
                        logLik = -BMest$n2llh / 2,
                        p = 1 + d,
                        aic = 2 * (1 + d) + BMest$n2llh,
                        vcov = BMest$vcov)
    } else {
        if(sum(lower > start) + sum(upper < start) > 0){
            stop("The starting value is not within the bounds of the parameter.")
        }

        minus2llh_sinvar <- function(logvalue) {
            if(model %in% c("BM", "trend")){
                prm[[1]] <- exp(logvalue)
            } else {
                prm[[2]] <- exp(logvalue)
            }
            loglik(prm, y, X)$n2llh
        }

        minus2llh <- function(logvalue, y) {
            if(model == "EB"){
                prm[[1]] <- logvalue[1]
            } else {
                prm[[1]] <- exp(logvalue[1])
            }
            if(measurement_error && !(model %in% c("BM", "trend"))){
                prm[[2]] <- exp(logvalue[2])
            }
            loglik(prm, y, X)$n2llh
        }

        if(lower[1] == upper[1] && !measurement_error){
            prm[[1]] <- lower[1]
            BMest <- loglik(prm, y, X)
        } else {
            if(model != "EB"){
                logstart <- log(start)
                loglower <- log(lower)
                logupper <- log(upper)
            } else {
                logstart <- start
                loglower <- lower
                logupper <- upper
                if(measurement_error){
                    logstart[2] <- log(start[2])
                    loglower[2] <- log(lower[2])
                    logupper[2] <- log(upper[2])
                }
            }

            if(!(model %in% c("BM", "trend")) && lower[1] != upper[1]){
                opt <- optim(logstart, fn = minus2llh, method = "L-BFGS-B",
                             lower = loglower, upper = logupper, y = y, ...)
                if(model == "EB"){
                    MLEvalue <- as.numeric(opt$par[1])
                } else {
                    MLEvalue <- as.numeric(exp(opt$par[1]))
                }
                prm[[1]] <- MLEvalue
                if(measurement_error){
                    MLEsigma2_error <- as.numeric(exp(opt$par[2]))
                    prm[[2]] <- MLEsigma2_error
                }
            } else {
                if(!(model %in% c("BM", "trend"))){
                    prm[[1]] <- lower[1]
                    logstart <- logstart[2]
                    loglower <- loglower[2]
                    logupper <- logupper[2]
                }
                opt <- optim(logstart, fn = minus2llh_sinvar, method = "L-BFGS-B",
                             lower = loglower, upper = logupper, ...)
                MLEsigma2_error <- as.numeric(exp(opt$par[1]))
                if(model %in% c("BM", "trend")){
                    prm[[1]] <- MLEsigma2_error
                } else {
                    prm[[2]] <- MLEsigma2_error
                }
            }

            matchbound <- NULL
            if(isTRUE(all.equal(prm[[1]], lower[1], tol = tol)) ||
               isTRUE(all.equal(prm[[1]], upper[1], tol = tol))){
                matchbound <- c(matchbound, 1)
                if((model %in% c("lambda", "kappa")) && (prm[[1]] == 1)) matchbound <- NULL
                if((model == "EB") && (prm[[1]] == 0)) matchbound <- NULL
            }
            if(length(lower) > 1 &&
               (isTRUE(all.equal(prm[[2]], lower[2], tol = tol)) ||
                isTRUE(all.equal(prm[[2]], upper[2], tol = tol)))){
                matchbound <- c(matchbound, 2)
            }
            if(length(matchbound) > 0){
                for(i in matchbound){
                    warning(paste("the estimation of", names(prm)[i],
                                  "matches the upper/lower bound for this parameter.\n                          You may change the bounds using options \"upper.bound\" and \"lower.bound\".\n"))
                }
            }
            BMest <- loglik(prm, y, X)
        }

        sigma2_errorhat <- 0
        if(measurement_error){
            sigma2_errorhat <- MLEsigma2_error * BMest$sigma2hat
        }
        if(model %in% OU){
            BMest$sigma2hat <- 2 * prm[[1]] * BMest$sigma2hat
        }
        results <- list(coefficients = BMest$betahat,
                        sigma2 = BMest$sigma2hat,
                        optpar = prm[[1]],
                        sigma2_error = sigma2_errorhat,
                        logLik = -BMest$n2llh / 2,
                        p = 2 + d,
                        aic = 2 * (2 + d) + BMest$n2llh,
                        vcov = BMest$vcov)
        if(model %in% c("BM", "trend")){
            results$optpar <- NULL
            results$p <- results$p - 1
            results$aic <- results$aic - 2
        }
        if(measurement_error){
            results$p <- results$p + 1
            results$aic <- results$aic + 2
        }
    }

    names(results$coefficients) <- colnames(X)
    colnames(results$vcov) <- colnames(X)
    rownames(results$vcov) <- colnames(X)
    results$fitted.values <- drop(X %*% results$coefficients)
    results$residuals <- y - results$fitted.values
    results$mean.tip.height <- Tmax
    results$y <- y
    results$X <- X
    results$n <- n
    results$d <- d
    results$formula <- formula
    results$call <- match.call()
    results$model <- model
    results$boot <- boot
    results$REML <- REML
    results$model.frame <- mf
    class(results) <- "phylolm"
    return(results)
}
