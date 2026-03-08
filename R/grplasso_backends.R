resolve_grplasso_backend <- function(opt) {

    backend <- opt$grplasso.backend
    if (is.null(backend) || is.na(backend[[1]]) || !nzchar(backend[[1]])) {
        return("cpp")
    }

    backend <- as.character(backend[[1]])
    if (!(backend %in% c("cpp", "vendor-r", "package"))) {
        stop("grplasso backend must be one of 'cpp', 'vendor-r', or 'package'.")
    }

    backend
}

grplasso_path_result <- function(coefficients, lambda, backend, converged = rep(TRUE, length(lambda))) {

    list(
        coefficients = as.matrix(coefficients),
        lambda = as.numeric(lambda),
        converged = as.logical(converged),
        call = paste0("grplasso_", backend)
    )
}

run_package_grplasso_lambdamax <- function(grpX, grpY, grpIdx) {

    suppressMessages(
        grplasso::lambdamax(
            as.matrix(grpX),
            as.numeric(grpY),
            model = grplasso::LinReg(),
            index = grpIdx,
            standardize = FALSE,
            center = FALSE
        )
    )
}

linreg_group_lasso_lambda_max_r <- function(grpX, grpY, grpIdx) {

    groups <- split(seq_along(grpIdx), grpIdx)
    grad0 <- drop(-2 * crossprod(as.matrix(grpX), as.numeric(grpY)))
    max(vapply(groups, function(ind) {
        sqrt(sum(grad0[ind]^2)) / sqrt(length(ind))
    }, numeric(1)))
}

linreg_group_lasso_lambda_max <- function(grpX, grpY, grpIdx, backend = c("cpp", "vendor-r", "package")) {

    backend <- match.arg(backend)

    if (backend == "package") {
        return(run_package_grplasso_lambdamax(grpX, grpY, grpIdx))
    }

    if (backend == "cpp") {
        out <- tryCatch(
            linreg_group_lasso_lambda_max_cpp(as.matrix(grpX), as.numeric(grpY), as.integer(grpIdx)),
            error = function(e) NULL
        )
        if (!is.null(out)) {
            return(as.numeric(out))
        }
    }

    linreg_group_lasso_lambda_max_r(grpX, grpY, grpIdx)
}

run_package_grplasso_path <- function(grpX, grpY, grpIdx, lambda, tol) {

    sol <- grplasso::grplasso(
        as.matrix(grpX),
        y = as.numeric(grpY),
        standardize = FALSE,
        center = FALSE,
        lambda = lambda,
        model = grplasso::LinReg(),
        index = grpIdx,
        control = grplasso::grpl.control(tol = tol, trace = 0, save.y = FALSE)
    )

    grplasso_path_result(sol$coefficients, lambda, "package", converged = sol$converged)
}

run_vendored_grplasso_path <- function(grpX, grpY, grpIdx, lambda, tol) {

    grpX <- as.matrix(grpX)
    grpY <- as.numeric(grpY)
    grpIdx <- as.integer(grpIdx)
    groups <- split(seq_along(grpIdx), grpIdx)
    x.groups <- lapply(groups, function(ind) grpX[, ind, drop = FALSE])
    sqrt.group.size <- sqrt(vapply(groups, length, integer(1)))
    nH.pen <- vapply(x.groups, function(xj) {
        max(1e-2, 2 * max(colSums(xj * xj)))
    }, numeric(1))

    nGroups <- length(groups)
    nLambda <- length(lambda)
    nCols <- ncol(grpX)
    beta.ls <- 0.5
    sigma.ls <- 0.1
    inner.loops <- 10L
    max.iter <- 500L
    sqrt.tol <- sqrt(tol)

    coef <- numeric(nCols)
    norms.pen <- numeric(nGroups)
    coef.path <- matrix(0, nrow = nCols, ncol = nLambda)
    converged <- rep(TRUE, nLambda)
    residual <- grpY
    rss <- sum(residual^2)

    for (pos in seq_along(lambda)) {
        lmbd <- lambda[[pos]]
        fn.val <- rss + sum(lmbd * sqrt.group.size * norms.pen)
        do.all <- FALSE
        d.fn <- 1
        d.par <- 1
        counter <- 1L
        iter.count <- 0L

        while (d.fn > tol || d.par > sqrt.tol || !do.all) {
            if (iter.count >= max.iter) {
                converged[[pos]] <- FALSE
                break
            }

            fn.val.old <- fn.val
            coef.old <- coef

            if (counter == 0L || counter > inner.loops) {
                do.all <- TRUE
                guessed.active <- seq_len(nGroups)
                counter <- 1L
            } else {
                guessed.active <- which(norms.pen != 0)
                if (length(guessed.active) == 0L) {
                    do.all <- TRUE
                    guessed.active <- seq_len(nGroups)
                } else {
                    do.all <- FALSE
                    counter <- counter + 1L
                }
            }

            if (do.all) {
                iter.count <- iter.count + 1L
            }

            start.pen <- rep(1, nGroups)

            for (j in guessed.active) {
                ind <- groups[[j]]
                xj <- x.groups[[j]]
                coef.ind <- coef[ind]
                ngrad <- drop(-2 * crossprod(xj, residual))
                nH <- nH.pen[[j]]
                cond <- -ngrad + nH * coef.ind
                cond.norm <- sqrt(sum(cond^2))
                border <- sqrt.group.size[[j]] * lmbd

                if (cond.norm > border) {
                    coef.target <- ((1 - border / cond.norm) / nH) * cond
                    d <- coef.target - coef.ind
                } else {
                    d <- -coef.ind
                }

                if (!any(d != 0)) {
                    norms.pen[[j]] <- sqrt(sum(coef.ind^2))
                    next
                }

                scale <- min(start.pen[[j]] / beta.ls, 1)
                xjd <- drop(xj %*% d)
                dot.res.xjd <- sum(residual * xjd)
                dot.xjd.xjd <- sum(xjd * xjd)
                coef.norm.old <- norms.pen[[j]]
                qh <- sum(ngrad * d) + border * (sqrt(sum((coef.ind + d)^2)) - coef.norm.old)
                rss.test <- rss - 2 * scale * dot.res.xjd + scale * scale * dot.xjd.xjd
                coef.test <- coef.ind + scale * d
                left <- rss.test + border * sqrt(sum(coef.test^2))
                right <- rss + border * coef.norm.old + sigma.ls * scale * qh

                while (left > right && scale > 1e-30) {
                    scale <- scale * beta.ls
                    rss.test <- rss - 2 * scale * dot.res.xjd + scale * scale * dot.xjd.xjd
                    coef.test <- coef.ind + scale * d
                    left <- rss.test + border * sqrt(sum(coef.test^2))
                    right <- rss + border * coef.norm.old + sigma.ls * scale * qh
                }

                if (scale <= 1e-30) {
                    start.pen[[j]] <- 1
                } else {
                    coef[ind] <- coef.test
                    residual <- residual - scale * xjd
                    rss <- rss.test
                    start.pen[[j]] <- scale
                }

                norms.pen[[j]] <- sqrt(sum(coef[ind]^2))
            }

            fn.val <- rss + sum(lmbd * sqrt.group.size * norms.pen)
            d.par <- max(abs(coef - coef.old) / (1 + abs(coef)))
            d.fn <- abs(fn.val.old - fn.val) / (1 + abs(fn.val))

            if (d.fn <= tol && d.par <= sqrt.tol) {
                counter <- 0L
            }
        }

        coef.path[, pos] <- coef
    }

    grplasso_path_result(coef.path, lambda, "vendor-r", converged = converged)
}

run_cpp_grplasso_path <- function(grpX, grpY, grpIdx, lambda, tol) {

    sol <- tryCatch(
        linreg_group_lasso_path_cpp(
            as.matrix(grpX),
            as.numeric(grpY),
            as.integer(grpIdx),
            as.numeric(lambda),
            tol,
            500L,
            10L,
            0.5,
            0.1,
            TRUE
        ),
        error = function(e) NULL
    )

    if (is.null(sol)) {
        return(run_vendored_grplasso_path(grpX, grpY, grpIdx, lambda, tol))
    }

    grplasso_path_result(sol$coefficients, lambda, "cpp", converged = sol$converged)
}

run_grplasso_path <- function(grpX, grpY, grpIdx, lambda, tol, backend = c("cpp", "vendor-r", "package")) {

    backend <- match.arg(backend)

    if (backend == "package") {
        return(run_package_grplasso_path(grpX, grpY, grpIdx, lambda, tol))
    }
    if (backend == "vendor-r") {
        return(run_vendored_grplasso_path(grpX, grpY, grpIdx, lambda, tol))
    }

    run_cpp_grplasso_path(grpX, grpY, grpIdx, lambda, tol)
}
