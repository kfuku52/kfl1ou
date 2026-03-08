# Vendored and adapted from the CRAN package phylolm 2.6.5 (GPL >= 2).
# Original package authors are credited in NOTICE.
# Modified for kfl1ou on 2026-03-08.

pruningwise.branching.times <- function(phy) {
    ## Calculates branching times = node ages, for an ultrametric tree
    ## in pruningwise order.
    xx <- numeric(phy$Nnode)
    nt <- length(phy$tip.label)
    interns <- which(phy$edge[, 2] > nt)
    for(i in rev(interns)){
        xx[phy$edge[i, 2] - nt] <- xx[phy$edge[i, 1] - nt] + phy$edge.length[i]
    }
    depth <- xx[phy$edge[1, 1] - nt] + phy$edge.length[1]
    xx <- depth - xx
    names(xx) <- if(is.null(phy$node.label)) {
        (nt + 1):(nt + phy$Nnode)
    } else {
        phy$node.label
    }
    return(xx)
}

pruningwise.distFromRoot <- function(phy) {
    ## Distance from root to all nodes, for a tree in pruningwise order.
    nt <- length(phy$tip.label)
    xx <- numeric(phy$Nnode + nt)
    for(i in length(phy$edge.length):1){
        xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
    }
    names(xx) <- if(is.null(phy$node.label)) {
        1:(nt + phy$Nnode)
    } else {
        c(phy$tip.label, phy$node.label)
    }
    return(xx)
}

transf.branch.lengths <- function(phy,
                                  model = c("BM", "OUrandomRoot", "OUfixedRoot",
                                            "lambda", "kappa", "delta", "EB",
                                            "trend"),
                                  parameters = NULL,
                                  check.pruningwise = TRUE,
                                  check.ultrametric = TRUE,
                                  D = NULL,
                                  check.names = TRUE) {
    if(!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
    model <- match.arg(model)
    if(model == "trend" && is.ultrametric(phy)){
        stop("the trend is unidentifiable for ultrametric trees.")
    }
    if(is.null(phy$edge.length)) stop("the tree has no branch lengths.")
    if(is.null(phy$tip.label)) stop("the tree has no tip labels.")
    if(check.pruningwise) phy <- reorder(phy, "pruningwise")
    n <- length(phy$tip.label)
    des <- phy$edge[, 2]
    externalEdge <- (des <= n)
    if(!is.null(phy$root.edge) && phy$root.edge > 0){
        stop("the tree is supposed to have no root edge (or of length 0).")
    }

    parameters.default <- c(0, 1, 1, 1, 0, 0)
    names(parameters.default) <- c("alpha", "lambda", "kappa", "delta",
                                   "rate", "sigma2_error")

    if(is.null(parameters)){
        parameters <- parameters.default
    } else if(!inherits(parameters, "list")) {
        stop("please specify parameters as a list().")
    } else {
        specified <- !c(is.null(parameters$alpha),
                        is.null(parameters$lambda),
                        is.null(parameters$kappa),
                        is.null(parameters$delta),
                        is.null(parameters$rate),
                        is.null(parameters$sigma2_error))
        parameters.user <- c(parameters$alpha,
                             parameters$lambda,
                             parameters$kappa,
                             parameters$delta,
                             parameters$rate,
                             parameters$sigma2_error)
        parameters <- parameters.default
        parameters[specified] <- parameters.user
    }

    p <- list(alpha = parameters[1],
              lambda = parameters[2],
              kappa = parameters[3],
              delta = parameters[4],
              rate = parameters[5],
              sigma2_error = parameters[6])

    root.edge <- 0
    diagWeight <- rep(1, n)
    errEdge <- rep(p$sigma2_error, n)

    if(model %in% c("BM", "trend")){
        edge.length <- phy$edge.length
    }

    OU <- c("OUrandomRoot", "OUfixedRoot")
    if(model %in% OU){
        if(check.ultrametric){
            D <- numeric(n)
            if(!is.ultrametric(phy)){
                dis <- pruningwise.distFromRoot(phy)
                D <- max(dis[1:n]) - dis[1:n]
                D <- D - mean(D)
                phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] +
                    D[des[externalEdge]]
            }
        } else {
            if(is.null(D)) stop("Provide D if you choose check.ultrametric=F")
            if(length(D) != n) stop("D should be a vector with one term for each tip in the tree")
            if(check.names){
                if(is.null(names(D))) stop("D is lacking names (tip labels)")
                ordr <- match(phy$tip.label, names(D))
                if(sum(is.na(ordr)) > 0) stop("names of D do not match the tree tip labels.")
                D <- D[ordr]
            }
        }
        times <- pruningwise.branching.times(phy)
        Tmax <- max(times)
        alpha <- p$alpha
        errEdge <- errEdge * exp(-2 * alpha * D[des[externalEdge]])
        if(model == "OUrandomRoot"){
            distFromRoot <- exp(-2 * alpha * times)
            d1 <- distFromRoot[phy$edge[, 1] - n]
            d2 <- numeric(nrow(phy$edge))
            d2[externalEdge] <- exp(-2 * alpha * D[des[externalEdge]])
            d2[!externalEdge] <- distFromRoot[des[!externalEdge] - n]
        }
        if(model == "OUfixedRoot"){
            distFromRoot <- exp(-2 * alpha * times) *
                (1 - exp(-2 * alpha * (Tmax - times)))
            d1 <- distFromRoot[phy$edge[, 1] - n]
            d2 <- numeric(nrow(phy$edge))
            d2[externalEdge] <- exp(-2 * alpha * D[des[externalEdge]]) *
                (1 - exp(-2 * alpha * (Tmax - D[des[externalEdge]])))
            d2[!externalEdge] <- distFromRoot[des[!externalEdge] - n]
        }
        edge.length <- d2 - d1
        root.edge <- min(distFromRoot)
        diagWeight <- exp(alpha * D)
    }

    if(model == "lambda"){
        lambda <- p$lambda
        distFromRoot <- pruningwise.distFromRoot(phy)
        edge.length <- phy$edge.length * lambda
        edge.length[externalEdge] <- edge.length[externalEdge] +
            (1 - lambda) * distFromRoot[des[externalEdge]]
    }

    if(model == "kappa"){
        edge.length <- phy$edge.length^p$kappa
    }

    if(model == "delta"){
        delta <- p$delta
        distFromRoot <- pruningwise.distFromRoot(phy)
        depth <- max(distFromRoot)
        edge.length <- (distFromRoot[des]^delta - distFromRoot[phy$edge[, 1]]^delta) *
            depth^(1 - delta)
    }

    if(model == "EB"){
        rate <- p$rate
        if(rate == 0){
            edge.length <- phy$edge.length
        } else {
            distFromRoot <- pruningwise.distFromRoot(phy)
            edge.length <- (exp(rate * distFromRoot[des]) -
                exp(rate * distFromRoot[phy$edge[, 1]])) / rate
        }
    }

    edge.length[externalEdge] <- edge.length[externalEdge] + errEdge
    phy$edge.length <- edge.length
    phy$root.edge <- root.edge
    names(diagWeight) <- phy$tip.label
    return(list(tree = phy, diagWeight = diagWeight))
}
