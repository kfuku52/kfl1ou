# 
#' Detects evolutionary shifts under an OU model
#'
#'This function takes in one or multiple traits, and automatically detects the phylogenetic placement and 
#'the magnitude of shifts in the evolution of these traits. The model assumes an Ornstein-Uhlenbeck process
#'whose parameters are estimated (adaptation `strength' \eqn{\alpha}{alpha} and drift variance \eqn{\sigma^2}{sigma^2}).
#'Instantaneous shifts in the optimal trait value affect the traits over time.
#'
#'@param tree ultrametric tree of class phylo with branch lengths, and edges in postorder.
#'@param Y trait vector/matrix without missing entries. The row names of the data must be in the same order as the tip labels.
#'@param max.nShifts upper bound for the number of shifts. The default value is half the number of tips.
#'@param criterion information criterion for model selection (see Details in \code{\link{configuration_ic}}).
#'@param root.model ancestral state model at the root.
#'@param candid.edges a vector of indices of candidate edges where the shifts may occur. If provided, shifts will only be allowed on these edges; otherwise all edges will be considered.
#'@param quietly logical. If FALSE, a basic summary of the progress and results is printed.
#'@param alpha.starting.value optional starting value for the optimization of the phylogenetic adaptation rate. 
#'@param alpha.upper optional upper bound for the phylogenetic adaptation rate. The default value is log(2) over the minimum branch length connected to tips. 
#'@param alpha.lower optional lower bound for the phylogenetic adaptation rate.
#'@param lars.alg sparse-path algorithm used in the univariate case. The
#' option name is kept for backward compatibility.
#'@param nCores maximum total CPU budget for \code{kfl1ou}. If \code{nCores=1}
#' then it runs sequentially. Otherwise, when fork-based parallelism is
#' available, \code{kfl1ou} may use up to \code{nCores} forked workers via
#' \code{mclapply}, while BLAS/OpenMP threads are limited when supported so the
#' overall computation respects this budget when possible.
#'@param rescale logical. If TRUE, the columns of the trait matrix are first rescaled so that all have the same l2-norm. If TRUE, the scores will be based on the rescale one.
#'@param edge.length.threshold minimum edge length that is considered non-zero. Branches with shorter length are considered as soft polytomies, disallowing shifts on such branches.
#'@param grp.delta internal (used when the data contain multiple traits). The input lambda sequence for the group lasso, in `grplasso', will be lambda.max*(0.5^seq(0, grp.seq.ub, grp.delta) ).
#'@param grp.seq.ub (used for multiple traits). The input lambda sequence for grplasso will be lambda.max*(0.5^seq(0, grp.seq.ub, grp.delta) ).
#'@param l1ou.options if provided, all the default values will be ignored. 
#'@return 
#' \item{Y}{input trait vector/matrix.}
#' \item{tree}{input tree.}
#' \item{shift.configuration}{estimated shift positions, i.e. vector of indices of edges where the estimated shifts occur.}
#' \item{shift.values}{estimates of the shift values.}
#' \item{shift.means}{estimates change of the expectation of the shift values}
#' \item{nShifts}{estimated number of shifts.}
#' \item{optima}{optimum values of the trait at tips. If the data are multivariate, this is a matrix where each row corresponds to a tip.}
#' \item{edge.optima}{optimum values of the trait on the edges. If the data are multivariate, this is a matrix where each row corresponds to an edge.}
#' \item{alpha}{maximum likelihood estimate(s) of the adaptation rate \eqn{\alpha}{alpha}, one per trait.}
#' \item{sigma2}{maximum likelihood estimate(s) of the variance rate \eqn{\sigma^2}{sigma^2}, one per trait.}
#' \item{mu}{fitted values, i.e. estimated trait means.}
#' \item{residuals}{residuals. These residuals are phylogenetically correlated.}
#' \item{score}{information criterion value of the estimated shift configuration.}
#' \item{profile}{list of shift configurations sorted by their ic scores.}
#' \item{l1ou.options}{list of options that were used.}
#'
#'@details
#'For information criteria: see \code{\link{configuration_ic}}.
#'
#'\code{nCores} is treated as the total CPU budget for \code{kfl1ou}, not only
#'as the number of forked worker processes. When \code{nCores > 1} and
#'fork-based parallelism is available, \code{kfl1ou} may use up to
#'\code{nCores} workers via \code{mclapply}. In that case BLAS/OpenMP threads
#'are reduced to 1 per process when supported, nested \code{kfl1ou} calls are run
#'sequentially, and the previous thread settings are restored on exit. On
#'platforms without fork-based parallelism, \code{kfl1ou} falls back to
#'sequential execution. For backward compatibility, fitted objects still use
#'class \code{"l1ou"} and store their options in the \code{l1ou.options}
#'component.
#'@examples
#' 
#' data(lizard.tree, lizard.traits)
#' keep <- lizard.tree$tip.label[1:15]
#' tree <- drop.tip(lizard.tree, setdiff(lizard.tree$tip.label, keep))
#' tree <- reorder(tree, "postorder")
#' Y <- lizard.traits[keep, 1]
#' lizard <- adjust_data(tree, Y)
#' eModel <- estimate_shift_configuration(
#'   lizard$tree, lizard$Y, criterion="AICc", max.nShifts=2
#' )
#' eModel
#'  
#' ## use up to two cores when available
#' eModel.par <- estimate_shift_configuration(
#'   lizard$tree, lizard$Y, criterion="AICc", max.nShifts=2, nCores=2
#' )
#' eModel.par$nShifts
#'
#' nEdges <- Nedge(lizard$tree) # total number of edges
#' ew <- rep(1, nEdges)  # to set default edge width of 1
#' if (length(eModel$shift.configuration) > 0) {
#'   ew[eModel$shift.configuration] <- 3   # widen edges with a shift
#' }
#' plot(eModel, cex=0.5, label.offset=0.02, edge.width=ew)
#'
#' # example to constrain the set of candidate branches with a shift
#' ce <- seq_len(min(5, Nedge(lizard$tree)))
#' eModel.ce <- estimate_shift_configuration(
#'   lizard$tree, lizard$Y, criterion="AICc", max.nShifts=2, candid.edges=ce
#' )
#' plot(eModel.ce, edge.ann.cex=0.7, cex=0.5, label.offset=0.02)
#'
#'@references
#'Mohammad Khabbazian, Ricardo Kriebel, Karl Rohe, and Cécile Ané (2016).
#' "Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models".
#' Methods in Ecology and Evolution. doi:10.1111/2041-210X.12534
#'
#'@export
estimate_shift_configuration <- function(tree, Y, 
           max.nShifts            = floor(length(tree$tip.label)/2), 
           criterion              = c("pBIC", "pBICess", "mBIC", "BIC", "AICc"), 
           root.model             = c("OUfixedRoot", "OUrandomRoot"),
           candid.edges           = NA,
           quietly                = TRUE,
           alpha.starting.value   = NA, 
           alpha.upper            = alpha_upper_bound(tree), 
           alpha.lower            = NA,
           lars.alg               = c("lasso", "stepwise"),
           nCores                 = 1,
           rescale                = TRUE,
           edge.length.threshold  = .Machine$double.eps,
           grp.delta              = 1/16,
           grp.seq.ub             = 5,
           measurement_error      = FALSE,
           input_error            = NULL,
           l1ou.options           = NA
     ){

    if (!inherits(tree, "phylo"))  stop("object \"tree\" is not of class \"phylo\".")
    if (is.null(tree$edge.length)) stop("the tree has no branch lengths.")
    if (is.null(tree$tip.label))   stop("the tree has no tip labels.")	
    if(!is.ultrametric(tree))      stop("the input tree is not ultrametric.")

    if( !identical(tree$edge, reorder(tree, "postorder")$edge))
        stop("the tree is not in postorder, use adjust_data function to reorder the tree!")

    Y <- as.matrix(Y)

    multivariate.missing <- FALSE
    if( any(is.na(Y)) ){
        if( ncol(Y) == 1){
            stop("some of the entries of the trait vector (Y) are missing.
                 you may use drop.tip to drop corresponding tips from the tree.\n")
        }else{
            ## the trait matrix can have some missing values as long as
            ## all the tips have at least one value.
            warning("some of the entries of the trait matrix (Y) are
                    missing.\n", immediate.=TRUE)
            if( length( which(rowSums(!is.na(Y))==0) ) != 0 ){
                stop("the trait matrix has a row with all missing values. you
                     may use drop.tip to drop corresponding tips from the tree.\n")
            }
            multivariate.missing <- TRUE
        }
    }

    if (!inherits(Y, "matrix")){
        Y <- as.matrix(Y)
    }

    if( nrow(Y) != length(tree$tip.label)){
       stop("the number of entries/rows of the trait vector/matrix (Y) 
            doesn't match the number of tips.\n") 
    }

    if( is.null(rownames(Y)) ){
        warning("no names provided for the trait(s) entries/rows. So it is assumed that 
                entries/rows match with the tip labels in the same order.\n", immediate.=TRUE)
        rownames(Y)  <- tree$tip.label
    } else{
        if( any(is.na(rownames(Y))) ){
            stop("some of the names in either trait vector/matrix or tree tip 
                 labels are unavailable.\n")
        }
    }

    if(!identical(rownames(Y), tree$tip.label)){
        diffres = setdiff(rownames(Y), tree$tip.label)
        if( length(diffres) > 0 ){
            cat(diffres)
            stop(" do(es) not exist in the tip labels of the input tree.\n")
        }
        diffres = setdiff(tree$tip.label, rownames(Y))
        if( length(diffres) > 0 ){
            cat(diffres)
            stop(" do(es) not exist in the input trait. You may want to use drop.tip(tree, setdiff(tree$tip.label,rownames(Y))) 
                 to drop extra tips in the tree.\n")
        }

        stop("the order of entries/rows of the trait vector/matrix (Y) does not matche 
             the order of the tip labels. Use adjust_data function to fix that.\n")
    }

    stopifnot(all(rownames(Y) == tree$tip.label))
    stopifnot(identical(rownames(Y), tree$tip.label))

    if( max.nShifts > length(tree$tip.label)){
        warning("max.nShifts should be a positive number less than number of tips. I set it to number of tips.\n")
        max.nShifts  <-  length(tree$tip.label)
    }
    if( max.nShifts < 0){
        warning("max.nShifts should be a positive number less than number of tips. I set it to 0.\n")
        max.nShifts  <- 1 
    }

    alpha.bounds <- sanitize_alpha_bounds(alpha.lower, alpha.upper)
    alpha.lower <- alpha.bounds$lower
    alpha.upper <- alpha.bounds$upper

    stopifnot( nCores > 0 )
    parallel.computing <- FALSE
    if(nCores>1){
        if(!l1ou_supports_multicore()){
            warning("fork-based parallel execution is unavailable; running sequentially.", immediate=TRUE)
            nCores <- 1
        }else{
            parallel.computing <- TRUE
        }
    }

    if (all(is.na(l1ou.options))){
        l1ou.options                   <- list()
        l1ou.options$nCores             <- nCores
        l1ou.options$parallel.computing <- parallel.computing
        ##NOTE: The saving_score functions are unprotected. To avoid race,
        ## I simply disable them in parallel mode until later that I figure out how to fix it.
        if(parallel.computing)
            l1ou.options$use.saved.scores  <- FALSE
        else
            l1ou.options$use.saved.scores  <- TRUE

        l1ou.options$max.nShifts       <- max.nShifts
        l1ou.options$criterion         <- match.arg(criterion)
        l1ou.options$lars.alg          <- match.arg(lars.alg)
        l1ou.options$root.model        <- match.arg(root.model)
        l1ou.options$quietly           <- quietly

        l1ou.options$alpha.starting.value   <- alpha.starting.value
        l1ou.options$alpha.upper.bound      <- alpha.upper
        l1ou.options$alpha.lower.bound      <- alpha.lower
        l1ou.options$edge.length.threshold  <- edge.length.threshold
        l1ou.options$rescale                <- rescale

        l1ou.options$grp.seq.ub   <- grp.seq.ub
        l1ou.options$grp.delta    <- grp.delta
        l1ou.options$candid.edges <- candid.edges
        l1ou.options$Z            <- generate_design_matrix(tree, "simpX")
        l1ou.options$grplasso.backend <- "cpp"
        l1ou.options$measurement_error <- measurement_error
        l1ou.options$input_error       <- normalize_input_error(tree, Y, input_error)
        ## each tree in tree.list represents a trait where the tips corresponding
        ## to NA values in the trail have been dropped.
        l1ou.options$multivariate.missing <- multivariate.missing
        if(multivariate.missing)
            l1ou.options$tree.list <- gen_tree_array(tree, Y)
        else
            l1ou.options$tree.list <- NULL
    } else {
        if( is.null(l1ou.options$measurement_error) ){
            l1ou.options$measurement_error <- FALSE
        }
        l1ou.options$input_error <- normalize_input_error(tree, Y, l1ou.options$input_error)
    }
    check_input_error_support(l1ou.options$measurement_error, l1ou.options$input_error)
    l1ou.options <- initialize_design_cache(tree, l1ou.options)
    l1ou.options <- initialize_fast_phylolm_cache(tree, l1ou.options)

    thread.limit <- resolve_l1ou_thread_limit(l1ou.options$nCores,
                                              l1ou.options$parallel.computing)
    return(with_l1ou_thread_limit(thread.limit, {
        if(length(l1ou.options$candid.edges)==0 || l1ou.options$max.nShifts == 0){
            l1ou.options$use.saved.scores <- FALSE
            return( fit_OU(tree, Y, shift.configuration=c(), l1ou.options=l1ou.options) )
        }

        if (l1ou.options$use.saved.scores) { erase_configuration_score_db() }

        if(!isTRUE(l1ou.options$quietly))
            cat("Starting first LASSO (alpha=0) to find a list of candidate configurations.\n")
        if (ncol(Y) == 1) {
            eModel1 = estimate_shift_configuration_known_alpha(tree, 
                Y, est.alpha = TRUE, opt = l1ou.options)
            if(!isTRUE(l1ou.options$quietly))
                cat("Starting second LASSO (alpha=",round(eModel1$alpha,2),") for another list of candidates.\n")
            eModel = estimate_shift_configuration_known_alpha(tree, 
                Y, alpha = eModel1$alpha, opt = l1ou.options)
            if (eModel$score > eModel1$score) 
                eModel = eModel1
        }
        if (ncol(Y) > 1) {
            eModel1 = estimate_shift_configuration_known_alpha_multivariate(tree, 
                Y, est.alpha = TRUE, opt = l1ou.options)
            if(!isTRUE(l1ou.options$quietly))
                cat("Starting second LASSO (alpha=",round(eModel1$alpha,2),") for another list of candidates.\n")
            candidate.configurations <- c(
                eModel1$.candidate.configurations,
                estimate_shift_configuration_known_alpha_multivariate(
                    tree, Y,
                    alpha = eModel1$alpha,
                    opt = l1ou.options,
                    candidate.only = TRUE
                )
            )
            eModel = estimate_shift_configuration_from_candidates(
                tree, Y,
                candidate.configurations = candidate.configurations,
                opt = l1ou.options
            )
            if (eModel$score > eModel1$score) 
                eModel = eModel1
        }

        eModel$profile = list_investigated_configs() 
        eModel$.candidate.configurations <- NULL

        if (l1ou.options$use.saved.scores) { erase_configuration_score_db() }
        if ( ! all( eModel$alpha < (alpha.upper - .Machine$double.eps)  ) )
		        warning('estimated alpha is too close to its upper bound. You may want to increase alpha.upper.\n')

        return(eModel)
    }))
}


estimate_shift_configuration_known_alpha <- function(tree, Y, alpha=0, est.alpha=FALSE, opt){  

    stopifnot( alpha >=0 )
    input.error <- get_trait_input_error(opt, 1, tree=tree)

    if ( est.alpha ){ ## BM model
        X   = cached_weighted_design_matrix(opt$Z, opt$edge.age, type="apprX")
        whitening.fit <- estimate_whitening_fit(tree, Y, alpha=0, est.alpha=TRUE,
                                                opt=opt, input_error=input.error)
        Cinvh   = t( sqrt_OU_covariance(tree, alpha=0, root.model = "OUfixedRoot",
                                        sigma2 = whitening.fit$sigma2,
                                        sigma2_error = whitening.fit$sigma2_error,
                                        input_error = input.error,
                                        check.order=F, check.ultrametric=F)$sqrtInvSigma ) 
    } else{           ## OU model
        X   = cached_weighted_design_matrix(opt$Z, opt$edge.age, type="orgX", alpha=alpha)
        whitening.fit <- estimate_whitening_fit(tree, Y, alpha=alpha, est.alpha=FALSE,
                                                opt=opt, input_error=input.error)
        Cinvh   = t( sqrt_OU_covariance(tree, alpha=alpha, root.model = opt$root.model,
                                        sigma2 = whitening.fit$sigma2,
                                        sigma2_error = whitening.fit$sigma2_error,
                                        input_error = input.error,
                                        check.order=F, check.ultrametric=F)$sqrtInvSigma ) 
    }

    if(!all(is.na(opt$candid.edges))){
        to.be.removed  = setdiff(1:length(tree$edge.length), opt$candid.edges)
    }else{
        to.be.removed  = c(length(tree$edge.length), which(tree$edge.length < opt$edge.length.threshold))
    }

    # NOTE: below: noise whitening, using contrasts: independent, same variance.
    # all but the last contrast have mean 0: Cinvh %*% column of ones = column of zeros, except last
    # last contrast: is about intercept, or ancestral state, should not be penalized.
    #YY  = as.matrix(Cinvh%*%Y)
    YY  = (Cinvh%*%Y)
    YY  = YY[-nrow(YY), ] 

    XX  = Cinvh%*%X

    nP  = ncol(XX)
    XX  = as.matrix(XX[,-to.be.removed])
    # NOTE: refer to the above note about whitening.
    XX  = XX[-nrow(XX), ]

    capture.output(
        sol.path <- run_univariate_sparse_path(XX, YY, opt)
    )

    Tmp = matrix(0, nrow(sol.path$beta), nP)
    Tmp[,-to.be.removed] = sol.path$beta
    sol.path$beta = Tmp

    result  = select_best_solution(tree, Y, sol.path, opt)
    eModel  = fit_OU_model(tree, Y, result$shift.configuration, opt)
    eModel$.candidate.configurations <- result$candidate.configurations

    if(!opt$quietly){
        print(eModel)
        cat("-------\n")
    }
    return(eModel)
}


estimate_shift_configuration_known_alpha_multivariate <- function(tree, Y, alpha=0, est.alpha=FALSE, opt,
                                                                  candidate.only=FALSE){

    stopifnot( alpha >=0 )

    if ( est.alpha == FALSE ){
        stopifnot(ncol(Y) == length(alpha))
    }

    Ymv <- Y
    input.error.mv <- opt$input_error
    if(opt$rescale==TRUE){
        rescaled <- rescale_matrix_and_error(Ymv, input.error.mv)
        Ymv <- rescaled$Y
        input.error.mv <- rescaled$input_error
    }

    nVariables    = ncol(Ymv)
    nEdges        = Nedge(tree)
    ##X             = generate_design_matrix(tree, "apprX")
    ##X             = cbind(X,1)
    ##ncolX         = ncol(X)
    ##to.be.removed = c(ncolX-1, which(tree$edge.length < opt$edge.length.threshold))

    if(!all(is.na(opt$candid.edges))){
        to.be.removed  = setdiff(1:length(tree$edge.length), opt$candid.edges)
    }else{
        to.be.removed = c(nEdges, which(tree$edge.length < opt$edge.length.threshold))
    }

    offset        = rep(nEdges*(0:(ncol(Ymv)-1)), each=length(to.be.removed))
    to.be.removed = rep(to.be.removed, ncol(Ymv)) + offset

    YY  = Ymv
    grpX.blocks = vector("list", ncol(Ymv))
    apprX <- NULL
    if(est.alpha){
        apprX <- cached_weighted_design_matrix(opt$Z, opt$edge.age, type="apprX")
    }
    for( i in 1:ncol(Ymv)){
        trait.input.error <- get_trait_input_error(opt, i, tree=tree, input_error=input.error.mv)

        if ( est.alpha == TRUE ){
            X   = apprX
            whitening.fit <- estimate_whitening_fit(tree, Ymv[,i,drop=FALSE], alpha=0,
                                                    est.alpha=TRUE, opt=opt,
                                                    input_error=trait.input.error)
            RE  = sqrt_OU_covariance(tree, root.model = "OUfixedRoot",
                                     alpha = 0,
                                     sigma2 = whitening.fit$sigma2,
                                     sigma2_error = whitening.fit$sigma2_error,
                                     input_error = trait.input.error,
                                     check.order=F, check.ultrametric=F)
        } else {
            X   = cached_weighted_design_matrix(opt$Z, opt$edge.age,
                                                type="orgX", alpha=alpha[[i]])
            whitening.fit <- estimate_whitening_fit(tree, Ymv[,i,drop=FALSE], alpha=alpha[[i]],
                                                    est.alpha=FALSE, opt=opt,
                                                    input_error=trait.input.error)
            RE  = sqrt_OU_covariance(tree,  root.model = opt$root.model,   
                                     alpha = alpha[[i]],
                                     sigma2 = whitening.fit$sigma2,
                                     sigma2_error = whitening.fit$sigma2_error,
                                     input_error = trait.input.error,
                                     check.order=F, check.ultrametric=F)
        }
        Cinvh   = t(RE$sqrtInvSigma) #\Sigma^{-1/2}

        y.ava        = !is.na(Ymv[,i])
        YY[y.ava, i] = as.matrix(Cinvh[y.ava, y.ava] %*% Ymv[y.ava, i])


        ##X     = cbin(Cinvh%*%X,1)
        ##grpX  = adiag(grpX, X)  
	XX      = Cinvh%*%X
	#NOTE: refer to whitening note above.
	XX      = XX[-nrow(XX), ]

	grpX.blocks[[i]] = XX
    }
    grpX = block_diag_matrix(grpX.blocks)
    #NOTE: refer to the whitening note above.
    YY  = YY[-nrow(YY), ] 

    np     = ncol(grpX)
    grpY   = c(YY)
    grpX   = as.matrix(grpX[,-to.be.removed])
    grpIdx = rep(1:ncol(X), ncol(Ymv))[-to.be.removed]
    ##grpIdx = rep(c(1:(ncol(X)-1),NA), ncol(Ymv))[-to.be.removed]


    ###NOTE: handling NAs in grpY
    grpY.ava = !is.na(grpY)
    grpY     = grpY[grpY.ava]
    grpX     = grpX[grpY.ava, ]

    grpX.col.nZero.idx = which(colSums(abs(grpX))!=0)
    grpX.nCol          = ncol(grpX)
    grpX               = grpX[,grpX.col.nZero.idx]
    grpIdx             = grpIdx[grpX.col.nZero.idx]

    sol = run_grplasso(grpX, grpY, nVariables, grpIdx, opt)

    Tmp                      = matrix(0, grpX.nCol, ncol(sol$coefficients))
    Tmp[grpX.col.nZero.idx,] = matrix(sol$coefficients)
    sol$coefficients         = Tmp
    ###end NOTE:
    
    Tmp                  = matrix(0, np, ncol(sol$coefficients))
    Tmp[-to.be.removed,] = matrix(sol$coefficients)
    sol$coefficients     = Tmp

    ##removing the intercept results
    #sol$coefficients     = sol$coefficients[-ncol(grpX), ]

    candidate.configurations <- collect_candidate_configurations(tree, Y, sol, opt)
    if(isTRUE(candidate.only)){
        return(candidate.configurations)
    }

    ### use the original Y
    result  = evaluate_candidate_configurations(tree, Y, candidate.configurations, opt=opt)
    eModel  = fit_OU_model(tree, Y, result$shift.configuration, opt=opt)
    eModel$.candidate.configurations <- result$candidate.configurations

    if(!opt$quietly){
        print(eModel)
        print("-------")
    }
    return(eModel)
}

estimate_shift_configuration_from_candidates <- function(tree, Y, candidate.configurations, opt){

    if(length(candidate.configurations) == 0){
        return(fit_OU_model(tree, Y, integer(), opt))
    }

    keys <- vapply(candidate.configurations, function(sc) {
        paste(sort(correct_unidentifiability(tree, sc, opt)), collapse=" ")
    }, character(1))
    keep <- !duplicated(keys)
    candidate.configurations <- candidate.configurations[keep]
    candidate.configurations <- Filter(function(sc) length(sc) <= opt$max.nShifts,
                                       candidate.configurations)

    if(length(candidate.configurations) == 0){
        return(fit_OU_model(tree, Y, integer(), opt))
    }

    search_ith_config <- function(sc){
        do_backward_correction(tree, Y, sc, opt)
    }

    if(opt$parallel.computing){
        all.res <- l1ou_mclapply(rev(candidate.configurations),
                                 FUN=search_ith_config,
                                 mc.cores=opt$nCores)
        all.res <- rev(all.res)
    } else{
        all.res <- lapply(candidate.configurations, search_ith_config)
    }

    scores <- vapply(all.res, function(x) x$score, numeric(1))
    result <- all.res[[which.min(scores)]]
    eModel <- fit_OU_model(tree, Y, result$shift.configuration, opt)
    eModel$.candidate.configurations <- candidate.configurations
    eModel
}



generate_design_matrix <- function(tree, type="apprX", alpha){
    stopifnot( is.ultrametric(tree) )
    stopifnot( sum( 1:length(tree$tip.label) %in% tree$edge[,1]) == 0)

    nTips  = length(tree$tip.label)
    nEdges = Nedge(tree)

    children <- tree_children_edges(tree)
    root2tip <- tree_root_to_tip_edge_paths(tree, children=children)
    X  = matrix(0, nTips, nEdges)
    for(i in seq_len(nTips)){
        X[i, root2tip[[i]]] = 1
    }

    if(type == "simpX"){
        return(X)
    }

    depths <- tree_node_depths(tree, children=children)
    Tval <- max(depths[seq_len(nTips)])
    edge.values <- Tval - depths[tree$edge[, 1]]

    if(type == "orgX"){
        edge.values <- 1 - exp(-alpha*edge.values)
    }else if(type == "apprX"){
        edge.values <- edge.values
    }else
        stop("Undefined design matrix type")

    return(sweep(X, 2, edge.values, FUN="*"))
}

select_best_solution <- function(tree, Y, sol.path, opt){

    candidate.configurations <- collect_candidate_configurations(tree, Y, sol.path, opt)
    evaluate_candidate_configurations(tree, Y, candidate.configurations, opt)
}

collect_candidate_configurations <- function(tree, Y, sol.path, opt){

    nSols = get_num_solutions(sol.path)
    stopifnot( nSols > 0 )

    all.shifts = numeric()
    prev.shift.configuration = NA

    candid.idx <- 1
    shift.configuration.list <- list()
    for (idx in 1:nSols) {

        shift.configuration = get_configuration_in_sol_path(sol.path, idx, Y)
        shift.configuration = correct_unidentifiability(tree, shift.configuration, opt)

        if ( length(shift.configuration) > opt$max.nShifts ){break}
        if ( setequal(shift.configuration, prev.shift.configuration) ){next}
        prev.shift.configuration = shift.configuration

        ## sorting shifts based on their age in the solution path
        all.shifts  = c(all.shifts, shift.configuration)
        freq.shifts = rep(0, length(shift.configuration))
        count = 1
        for( s in shift.configuration){
            freq.shifts[[count]] = length( which(all.shifts == s) )
            count = count + 1
        }
        names(shift.configuration)  <- freq.shifts
        shift.configuration <- shift.configuration[order(names(shift.configuration), decreasing=TRUE)]
        shift.configuration.list[[candid.idx]] <- shift.configuration
        candid.idx <- candid.idx + 1
    }

    shift.configuration.list
}

evaluate_candidate_configurations <- function(tree, Y, candidate.configurations, opt){

    min.score = Inf
    best.shift.configuration <- numeric()

    if(length(candidate.configurations) == 0){
        return(list(score=min.score,
                    shift.configuration=best.shift.configuration,
                    candidate.configurations=candidate.configurations))
    }

    search_ith_config <- function(sc){
        res <- do_backward_correction(tree, Y, sc, opt)
        return(res)
    }

    if(opt$parallel.computing){
        all.res <- l1ou_mclapply(rev(candidate.configurations),
                                 FUN=search_ith_config,
                                 mc.cores=opt$nCores)
        for (i in length(all.res):1 ){
            res <- all.res[[i]] 
            if (min.score > res$score){
                min.score = res$score
                best.shift.configuration = res$shift.configuration
            }
        }
    }else{
        for (i in 1:length(candidate.configurations) ){
            res <- search_ith_config(candidate.configurations[[i]])
            if (min.score > res$score){
                min.score = res$score
                best.shift.configuration = res$shift.configuration
            }
        }
    }

    return ( list(score=min.score,
                  shift.configuration=best.shift.configuration,
                  candidate.configurations=candidate.configurations) )
}

do_backward_correction <- function(tree, Y, shift.configuration, opt){

    org.score = cmp_model_score(tree, Y, shift.configuration, opt)

    if( length(shift.configuration) < 3 ) { 
        return(list(score=org.score, shift.configuration=shift.configuration)) 
    }  

    #nShifts = length(shift.configuration)
    #removal.candids = shift.configuration[1:(nShifts-1)]
    #for( sp in removal.candids )
    for(sp in shift.configuration)
    {
        new.configuration = setdiff(shift.configuration, sp)
        new.score         = cmp_model_score(tree, Y, new.configuration, opt)      
        if ( new.score < org.score ){
            shift.configuration = new.configuration
            org.score           = new.score
        }
    }

    return(list(score=org.score, shift.configuration=shift.configuration))
}


#
#' Computes the information criterion score for a given configuration
#'
#'@param tree ultrametric tree of class phylo, with branch lengths, and edges in postorder.
#'@param Y trait vector/matrix without missing entries. The row names of the data must be in the same order as the tip labels.
#'@param shift.configuration shift positions, i.e. vector of indices of the edges where the shifts occur.
#'@param criterion an information criterion (see Details).
#'@param root.model an ancestral state model at the root.
#'@param alpha.starting.value optional starting value for the optimization of the phylogenetic adaptation rate. 
#'@param alpha.upper optional upper bound for the phylogenetic adaptation rate. The default value is log(2) over the minimum length of external branches, corresponding to a half life greater or equal to the minimum external branch length.
#'@param alpha.lower optional lower bound for the phylogenetic adaptation rate.
#'@param fit.OU.model logical. If TRUE, it returns an object of class l1ou with all the parameters estimated.
#'@param l1ou.options if provided, all the default values will be ignored. 
#'
#'@return Information criterion value of the given shift configuration.
#'
#'@details
#'AICc gives the usual small-sample size modification of AIC. 
#'BIC gives the usual Bayesian information criterion, here penalizing each shift as 2 parameters. 
#'mBIC is the modified BIC proposed by Ho and Ané (2014).
#'pBIC is the phylogenetic BIC for shifts proposed by Khabbazian et al.
#'pBICess is a version of pBIC where the determinant term is replaced by a sum of the log of effective sample sizes (ESS), similar to the ESS proposed by Ané (2008). 
#' 
#'@examples
#' 
#' data(lizard.tree, lizard.traits)
#' keep <- lizard.tree$tip.label[1:15]
#' tree <- drop.tip(lizard.tree, setdiff(lizard.tree$tip.label, keep))
#' tree <- reorder(tree, "postorder")
#' lizard <- adjust_data(tree, lizard.traits[keep, 1])
#' eModel <- estimate_shift_configuration(
#'   lizard$tree, lizard$Y, criterion="AICc", max.nShifts=2
#' )
#' configuration_ic(lizard$tree, eModel$Y, eModel$shift.configuration, criterion="pBIC")
#'
#' ### building l1ou object out of the second best score 
#' candidate <- eModel$profile$configurations[[min(2, length(eModel$profile$configurations))]]
#' eModel2 <- configuration_ic(
#'   eModel$tree, eModel$Y, candidate,
#'   fit.OU.model=TRUE, l1ou.options=eModel$l1ou.options
#' )
#' plot(eModel2)
#'
#'@seealso \code{\link{estimate_shift_configuration}} \code{\link{adjust_data}}
#'
#'@references
#'Cécile Ané, 2008. "Analysis of comparative data with hierarchical autocorrelation". Annals of Applied Statistics 2(3):1078-1102.
#'
#'Ho, L. S. T. and Ané, C. 2014.  "Intrinsic inference difficulties for trait evolution with Ornstein-Uhlenbeck models". Methods in Ecology and Evolution. 5(11):1133-1146.
#'
#'Mohammad Khabbazian, Ricardo Kriebel, Karl Rohe, and Cécile Ané (2016).
#' "Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models".
#' Methods in Ecology and Evolution. doi:10.1111/2041-210X.12534
#'
#'@export
configuration_ic <- function(tree, Y, shift.configuration, 
                     criterion    = c("pBIC", "pBICess", "mBIC", "BIC", "AICc"), 
                     root.model   = c("OUfixedRoot", "OUrandomRoot"),
                     alpha.starting.value = NA,
                     alpha.upper  = alpha_upper_bound(tree), 
                     alpha.lower  = NA,
                     measurement_error = FALSE,
                     input_error = NULL,
                     fit.OU.model = FALSE, 
                     l1ou.options = NA
                   ){

    if (!inherits(tree, "phylo"))  stop("object \"tree\" is not of class \"phylo\".")
    if( !identical(tree$edge, reorder(tree, "postorder")$edge))
        stop("the input phylogenetic tree is not in postorder. Use adjust_data function.")

    Y  = as.matrix(Y)
    if(!identical(rownames(Y), tree$tip.label)) stop("rownames of Y and tree$tip.label are not identical.")


    opt = list()
    if(!all(is.na(l1ou.options))){
        opt = l1ou.options
    }else{
        opt$criterion            <- match.arg(criterion)
        opt$root.model           <- match.arg(root.model)
        opt$quietly              <- TRUE
        opt$alpha.starting.value <- alpha.starting.value
        opt$alpha.upper.bound    <- alpha.upper
        opt$alpha.lower.bound    <- alpha.lower
        opt$Z                    <- generate_design_matrix(tree, "simpX")
        opt$multivariate.missing <- FALSE
        opt$use.saved.scores     <- FALSE
        opt$measurement_error    <- measurement_error
        opt$input_error          <- normalize_input_error(tree, Y, input_error)
    }
    if( is.null(opt$measurement_error) ){
        opt$measurement_error <- FALSE
    }
    if( is.null(opt$quietly) ){
        opt$quietly <- TRUE
    }
    alpha.bounds <- sanitize_alpha_bounds(opt$alpha.lower.bound, opt$alpha.upper.bound)
    opt$alpha.lower.bound <- alpha.bounds$lower
    opt$alpha.upper.bound <- alpha.bounds$upper
    opt$input_error <- normalize_input_error(tree, Y, opt$input_error)
    check_input_error_support(opt$measurement_error, opt$input_error)
    opt <- initialize_design_cache(tree, opt)
    opt <- initialize_fast_phylolm_cache(tree, opt)

    thread.limit <- resolve_l1ou_thread_limit(ifelse(is.null(opt$nCores), 1L, opt$nCores),
                                              FALSE)
    return(with_l1ou_thread_limit(thread.limit, {
        s.c = correct_unidentifiability(tree, shift.configuration, opt)
        if( length(s.c) != length(shift.configuration) )
            stop(paste0("the input shift configuration is not a parsimony configuration. 
                        For instance,\n", s.c, "\n is an alternative configuration with fewer shifts."))

        if( fit.OU.model ){
            ##TODO: convergent evolution.
            eModel = fit_OU_model(tree, Y, shift.configuration, opt)
            return(eModel)
        }
        score = cmp_model_score(tree, Y, shift.configuration, opt)
        return(score)
    }))
}


#
#' Fits an OU model based on a given configuration
#'
#'@param tree ultrametric tree of class phylo, with branch lengths, and edges in postorder.
#'@param Y trait vector/matrix without missing entries. The row names of the data must be in the same order as the tip labels.
#'@param shift.configuration shift positions, i.e. vector of indices of the edges where the shifts occur.
#'@param criterion an information criterion (see Details).
#'@param root.model model for the ancestral state at the root.
#'@param cr.regimes optional list of convergent regimes. Each element contains
#' the background regime (`0`) or a set of shift edges that should share the
#' same optimum.
#'@param alpha.starting.value optional starting value for the optimization of the phylogenetic adaptation rate. 
#'@param alpha.upper optional upper bound for the phylogenetic adaptation rate. The default value is log(2) over the minimum length of external branches, corresponding to a half life greater or equal to the minimum external branch length.
#'@param alpha.lower optional lower bound for the phylogenetic adaptation rate.
#'@param l1ou.options if provided, all the default values will be ignored. 
#'
#'@return an object of class l1ou similar to \code{\link{estimate_shift_configuration}}.
#'
#'@details
#'AICc gives the usual small-sample size modification of AIC. 
#'BIC gives the usual Bayesian information criterion, here penalizing each shift as 2 parameters. 
#'mBIC is the modified BIC proposed by Ho and Ané (2014).
#'pBIC is the phylogenetic BIC for shifts proposed by Khabbazian et al.
#'pBICess is a version of pBIC where the determinant term is replaced by a sum of the log of effective sample sizes (ESS), similar to the ESS proposed by Ané (2008). 
#' 
#'@examples
#' 
#' data(lizard.tree, lizard.traits)
#' keep <- lizard.tree$tip.label[1:15]
#' tree <- drop.tip(lizard.tree, setdiff(lizard.tree$tip.label, keep))
#' tree <- reorder(tree, "postorder")
#' lizard <- adjust_data(tree, lizard.traits[keep, 1])
#' eModel <- estimate_shift_configuration(
#'   lizard$tree, lizard$Y, criterion="AICc", max.nShifts=2
#' )
#'
#' ### building l1ou object out of the second best score 
#' candidate <- eModel$profile$configurations[[min(2, length(eModel$profile$configurations))]]
#' eModel2 <- fit_OU(
#'   eModel$tree, eModel$Y, candidate, l1ou.options=eModel$l1ou.options
#' )
#' plot(eModel2)
#' 
#' ### hypothesis testing
#'
#' shift.config <- 1L
#' hModel <- fit_OU(lizard$tree, lizard$Y, shift.config, criterion="AICc")
#' summary(hModel)
#'
#'@seealso \code{\link{estimate_shift_configuration}} \code{\link{adjust_data}}
#'
#'@references
#'Cécile Ané, 2008. "Analysis of comparative data with hierarchical autocorrelation". Annals of Applied Statistics 2(3):1078-1102.
#'
#'Ho, L. S. T. and Ané, C. 2014.  "Intrinsic inference difficulties for trait evolution with Ornstein-Uhlenbeck models". Methods in Ecology and Evolution. 5(11):1133-1146.
#'
#'Mohammad Khabbazian, Ricardo Kriebel, Karl Rohe, and Cécile Ané (2016).
#' "Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models".
#' Methods in Ecology and Evolution. doi:10.1111/2041-210X.12534
#'
#'@export
fit_OU <- function(tree, Y, shift.configuration, 
                     criterion    = c("pBIC", "pBICess", "mBIC", "BIC", "AICc"), 
                     root.model   = c("OUfixedRoot", "OUrandomRoot"),
                     cr.regimes   = NULL,
                     alpha.starting.value = NA,
                     alpha.upper  = alpha_upper_bound(tree), 
                     alpha.lower  = NA,
                     measurement_error = FALSE,
                     input_error = NULL,
                     l1ou.options = NA
                   ){

    if (!inherits(tree, "phylo"))  stop("object \"tree\" is not of class \"phylo\".")
    if( !identical(tree$edge, reorder(tree, "postorder")$edge))
        stop("the input phylogenetic tree is not in postorder. Use adjust_data function.")

    Y  = as.matrix(Y)
    if(!identical(rownames(Y), tree$tip.label)) stop("rownames of Y and tree$tip.label are not identical.")

    opt = list()
    if(!all(is.na(l1ou.options))){
        opt = l1ou.options
    }else{
        opt$criterion            <- match.arg(criterion)
        opt$root.model           <- match.arg(root.model)
        opt$quietly              <- TRUE
        opt$alpha.starting.value <- alpha.starting.value
        opt$alpha.upper.bound    <- alpha.upper
        opt$alpha.lower.bound    <- alpha.lower
        opt$Z                    <- generate_design_matrix(tree, "simpX")
        opt$multivariate.missing <- FALSE
        opt$use.saved.scores     <- FALSE
        opt$measurement_error    <- measurement_error
        opt$input_error          <- normalize_input_error(tree, Y, input_error)
    }
    if( is.null(opt$measurement_error) ){
        opt$measurement_error <- FALSE
    }
    if( is.null(opt$quietly) ){
        opt$quietly <- TRUE
    }
    alpha.bounds <- sanitize_alpha_bounds(opt$alpha.lower.bound, opt$alpha.upper.bound)
    opt$alpha.lower.bound <- alpha.bounds$lower
    opt$alpha.upper.bound <- alpha.bounds$upper
    opt$input_error <- normalize_input_error(tree, Y, opt$input_error)
    check_input_error_support(opt$measurement_error, opt$input_error)
    opt <- initialize_design_cache(tree, opt)
    opt <- initialize_fast_phylolm_cache(tree, opt)

    thread.limit <- resolve_l1ou_thread_limit(ifelse(is.null(opt$nCores), 1L, opt$nCores),
                                              FALSE)
    return(with_l1ou_thread_limit(thread.limit, {
        s.c = correct_unidentifiability(tree, shift.configuration, opt)
        if( length(s.c) != length(shift.configuration) )
            stop(paste0("the input shift configuration is not parsimonious. For instance, shifts on these edges:\n",
                        s.c, "provides an alternative equivalent configuration with fewer shifts."))

         eModel = fit_OU_model(tree, Y, shift.configuration, opt)
         if(!is.null(cr.regimes) ){
            if( isTRUE(opt$measurement_error) || !is.null(opt$input_error) ){
                stop("convergent regime fitting does not yet support measurement_error or input_error.")
            }

            if( !( 0 %in% unlist(cr.regimes) ) ){
                stop("background/intercept is not included in the regimes! Represent the background by \"0\".")
            }
            if( !(identical(sort(c(0,shift.configuration)), sort(unlist(cr.regimes)) ) ) ){
                stop("convergent regimes do not match with the shift positions.")
            }
            cr.score <- cmp_model_score_CR(tree, Y, regimes=cr.regimes,
                               alpha=eModel$alpha, opt=opt)
            eModel$cr.score <- cr.score
            for(idx in 1:length(cr.regimes)){
                if( 0 %in% cr.regimes[[idx]] ){
                    names(eModel$shift.configuration)[which(eModel$shift.configuration %in% cr.regimes[[idx]])] <- 0
                }else{
                    names(eModel$shift.configuration)[which(eModel$shift.configuration %in% cr.regimes[[idx]])] <- idx
                }
            }
         }
         return(eModel)
    }))
}

fit_OU_model <- function(tree, Y, shift.configuration, opt){

    Y      = as.matrix(Y)
    nEdges = Nedge(tree)
    nTips  = length(tree$tip.label)

    resi = mu = optima = matrix(data=NA, nrow=nTips, ncol=ncol(Y))
    if(length(shift.configuration) > 0){
        shift.values <- matrix(NA_real_, nrow=length(shift.configuration), ncol=ncol(Y))
        shift.means <- matrix(NA_real_, nrow=length(shift.configuration), ncol=ncol(Y))
    } else{
        shift.values <- numeric()
        shift.means <- numeric()
    }
    # edge.optima = matrix(NA, nEdges, ncol(Y))
    intercept = alpha = sigma2 = sigma2_error = rep(NA, ncol(Y))
    logLik = rep(NA, ncol(Y))

    for(i in 1:ncol(Y)){
        s.c <- c()
        if(!is.null(opt$tree.list)){
            tr <- opt$tree.list[[i]]
            prepared.tree <- opt$prepared.tree.list[[i]]
            edge.age <- tr$edge.age
            y  <- as.matrix(Y[!is.na(Y[,i]), i])
            input.error <- get_trait_input_error(opt, i, tree=tr, available=!is.na(Y[,i]))
            if(length(shift.configuration) > 0){
                augmented.s.c <- tr$old.order[shift.configuration] # index of shift edges in pruned tree, see gen_tree_array
                s.c <- as.integer(augmented.s.c[!is.na(augmented.s.c)])
            } else{
                augmented.s.c <- c()
            }
        } else{
            tr  <- tree
            prepared.tree <- opt$prepared.tree
            edge.age <- opt$edge.age
            y   <- as.matrix(Y[,i])
            s.c <- shift.configuration
            input.error <- get_trait_input_error(opt, i, tree=tr)
            augmented.s.c <- shift.configuration
        }

        nShifts = length(s.c)

        fit <- my_phylolm_interface(tr, y, s.c, opt, input_error=input.error,
                                    prepared.tree = prepared.tree)
        if ( all(is.na(fit)) ){
            stop("model score is NA in fit_OU_model function! 
		 This should not happen. Please set quietly to false to see the reason.")
        }

        alpha[i]  <- fit$optpar
        sigma2[i] <- fit$sigma2
        sigma2_error[i] <- fit$sigma2_error
        logLik[i] <- fit$logLik

        ## E[Y] and residuals: Y-EY
        mu[!is.na(Y[,i]), i]   = fit$fitted.values
        resi[!is.na(Y[,i]), i] = fit$residuals
        intercept[i] = fit$coefficients[[1]] # = y0 * e^-T + theta0_root * (1-e^-T), assumes ultrametric tree

        ## Now we have the alpha hat and we can form the true design matrix
        if( nShifts > 0 ){
            scale.values <- edge_scaling_from_cache(edge.age, type="orgX", alpha=alpha[i])[s.c]^-1
            fit$coefficients[2:(nShifts+1)] <- scale.values * fit$coefficients[2:(nShifts+1)]
        }

        if( length(shift.configuration) > 0 && nShifts > 0 ){
            visible.idx <- which(!is.na(augmented.s.c))
            shift.values[visible.idx, i] <- fit$coefficients[2:(nShifts+1)]
            shift.means[visible.idx, i] <- fit$coefficients[2:(nShifts+1)] / scale.values
        }

        optima.tmp = rep(fit$coefficients[[1]], nTips)  # optima at the tips for one trait
        if( length(shift.configuration) > 0 && nShifts > 0 ){
            visible.idx <- which(!is.na(augmented.s.c))
            if(length(visible.idx) > 0){
                optima.tmp <- optima.tmp + drop(
                    opt$Z[, shift.configuration[visible.idx], drop=FALSE] %*%
                        fit$coefficients[2:(nShifts+1)]
                )
            }
        }

        optima[,i] <- optima.tmp
    }

    rownames(optima) <- tree$tip.label

    ##NOTE: it reads the score from the database 
    ## and do not recompute the score. So it doesn't have any overhead.
    score = cmp_model_score (tree, Y, shift.configuration, opt) 
    model.opt <- opt
    model.opt$prepared.tree <- NULL
    model.opt$prepared.tree.list <- NULL

    model = list(
                 Y                   = Y, 
                 tree                = tree,
                 shift.configuration = shift.configuration, 
                 shift.values       = shift.values,
                 shift.means        = shift.means, 
                 nShifts             = length(shift.configuration), 
                 optima              = optima, 
                 alpha               = alpha, 
                 sigma2              = sigma2, 
                 sigma2_error        = sigma2_error,
                 intercept           = intercept, 
                 mu                  = mu, 
                 residuals           = resi,
                 score               = score,
                 logLik              = logLik,
                 l1ou.options        = model.opt) 

    class(model) <- "l1ou"
    return( model )
}

cmp_model_score <-function(tree, Y, shift.configuration, opt){

    ##TODO: optimize it.
    shift.configuration <- correct_unidentifiability(tree, shift.configuration, opt)

    if(opt$use.saved.scores){ ##if it's been already computed
        score <- get_configuration_score_from_list(shift.configuration)
        if(!is.na(score)){ return(score) }
    }

    Y  <- as.matrix(Y)
    ic <- opt$criterion

    res <- NA
    if( ic == "AIC"){
        stop("undefined")
    } else if( ic == "BIC"){
        res <- cmp_BIC(tree, Y, shift.configuration, opt)
    } else if( ic == "AICc"){
        res <- cmp_AICc(tree, Y, shift.configuration, opt)
    } else if( ic == "mBIC"){
        res <- cmp_mBIC(tree, Y, shift.configuration, opt) 
    } else if( ic == "pBICess"){
        res <- cmp_pBICess(tree, Y, shift.configuration, opt) 
    } else if(ic == "pBIC"){
        res <- cmp_pBIC(tree, Y, shift.configuration, opt) 
    } 

    if(all(is.na(res))){return(Inf)}

    score <- res$score
    if(opt$use.saved.scores){
        add_configuration_score_to_list(shift.configuration, score,
             paste0(c(res$sigma2/(2*res$alpha),res$logLik),collapse=" "))
    }
    return(score)
}

normalize_input_error <- function(tree, Y, input_error){

    if( is.null(input_error) ){
        return(NULL)
    }

    Y <- as.matrix(Y)
    if( inherits(input_error, "data.frame") ){
        input_error <- as.matrix(input_error)
    }

    if( is.null(dim(input_error)) ){
        if( length(input_error) != nrow(Y) ){
            stop("input_error must have one value per tip, or a matrix with one column per trait.")
        }
        err.mat <- matrix(input_error, nrow=nrow(Y), ncol=ncol(Y))
        if( !is.null(names(input_error)) ){
            rownames(err.mat) <- names(input_error)
        }
    } else {
        err.mat <- as.matrix(input_error)
        if( nrow(err.mat) != nrow(Y) ){
            stop("input_error must have the same number of rows as Y.")
        }
        if( ncol(err.mat) == 1 && ncol(Y) > 1 ){
            err.mat <- err.mat[, rep(1, ncol(Y)), drop=FALSE]
        }
        if( ncol(err.mat) != ncol(Y) ){
            stop("input_error must have one column per trait, or a single column shared across traits.")
        }
    }

    if( !is.null(rownames(err.mat)) ){
        o <- match(tree$tip.label, rownames(err.mat))
        if( any(is.na(o)) ){
            stop("input_error row names do not match tree tip labels.")
        }
        err.mat <- err.mat[o, , drop=FALSE]
    } else {
        rownames(err.mat) <- tree$tip.label
    }

    if( ncol(err.mat) > 1 && !all(is.null(colnames(Y))) ){
        if( !is.null(colnames(err.mat)) ){
            o <- match(colnames(Y), colnames(err.mat))
            if( any(is.na(o)) ){
                stop("input_error column names do not match the trait names.")
            }
            err.mat <- err.mat[, o, drop=FALSE]
        } else {
            colnames(err.mat) <- colnames(Y)
        }
    }

    if( any(err.mat < 0, na.rm=TRUE) ){
        stop("input_error must be non-negative.")
    }
    if( any(is.na(err.mat) & !is.na(Y)) ){
        stop("input_error has missing values for observed traits.")
    }

    return(err.mat)
}

sanitize_alpha_bounds <- function(alpha.lower, alpha.upper){

    lower <- ifelse(is.null(alpha.lower), NA_real_, as.numeric(alpha.lower[[1]]))
    upper <- ifelse(is.null(alpha.upper), NA_real_, as.numeric(alpha.upper[[1]]))

    if(!is.na(upper) && upper <= 0){
        stop("alpha.upper must be strictly positive.\n")
    }
    if(!is.na(lower) && lower < 0){
        warning("alpha.lower must be greater than zero. I set it to zero.\n")
        lower <- 0
    }
    if(!is.na(lower) && !is.na(upper) && upper < lower){
        warning("alpha.upper must be equal or greater than alpha.lower. I set them equal.\n")
        upper <- lower
    }

    return(list(lower=lower, upper=upper))
}

get_trait_input_error <- function(opt, idx, tree=NULL, available=NULL, input_error=opt$input_error){

    if( is.null(input_error) ){
        return(NULL)
    }

    col.idx <- ifelse(ncol(input_error) == 1, 1, idx)
    err <- input_error[, col.idx]
    err.names <- rownames(input_error)
    if( !is.null(available) ){
        err <- err[available]
        err.names <- err.names[available]
    }
    names(err) <- err.names

    if( !is.null(tree) ){
        o <- match(tree$tip.label, names(err))
        if( any(is.na(o)) ){
            stop("input_error names do not match the trait tree.")
        }
        err <- err[o]
        names(err) <- tree$tip.label
    }

    return(err)
}

prepare_trait_phylolm_data <- function(tree, Y, input_error=NULL){

    Y <- as.matrix(Y)
    y.names <- rownames(Y)
    if( is.null(y.names) ){
        y.names <- tree$tip.label
        rownames(Y) <- y.names
    }

    keep.data <- !is.na(Y[, 1])
    if( !all(keep.data) ){
        Y <- Y[keep.data, , drop=FALSE]
        y.names <- rownames(Y)
        if( !is.null(input_error) ){
            if( is.null(names(input_error)) ){
                stop("input_error must have names when matching observed traits to the tree.")
            }
            input_error <- input_error[y.names]
        }
    }

    if( !identical(y.names, tree$tip.label) || nrow(Y) != length(tree$tip.label) ){
        keep <- intersect(tree$tip.label, y.names)
        if( length(keep) < 2 ){
            stop("not enough taxa remain after matching the trait data to the tree.")
        }
        tree <- drop.tip(tree, setdiff(tree$tip.label, keep))
        tree <- reorder(tree, "postorder")
        Y <- Y[tree$tip.label, , drop=FALSE]
    }

    if( !is.null(input_error) ){
        if( is.null(names(input_error)) ){
            stop("input_error must have names when trait data are matched to a subset of tips.")
        }
        input_error <- input_error[tree$tip.label]
    }

    return(list(tree=tree, Y=Y, input_error=input_error))
}

estimate_whitening_fit <- function(tree, Y, alpha=0, est.alpha=FALSE, opt, input_error=NULL){

    if( !isTRUE(opt$measurement_error) && is.null(input_error) ){
        return(list(sigma2=1, sigma2_error=0))
    }

    prep <- prepare_trait_phylolm_data(tree, Y, input_error=input_error)
    yy <- as.numeric(prep$Y[, 1])
    names(yy) <- rownames(prep$Y)
    if( requires_phylolm_input_error_support(opt$measurement_error, prep$input_error) ){
        stop("input_error is not supported when measurement_error=TRUE in the current kfl1ou likelihood solver. Set measurement_error=FALSE.")
    }
    if( can_use_dense_input_error_fit(opt$measurement_error, prep$input_error) ){
        return(
            dense_known_input_error_gls_fit(
                prep$tree,
                prep$Y,
                preds = matrix(1, nrow=nrow(prep$Y), ncol=1),
                model = ifelse(est.alpha, "BM", opt$root.model),
                lower.bound = ifelse(est.alpha, NA_real_, alpha),
                upper.bound = ifelse(est.alpha, NA_real_, alpha),
                starting.value = ifelse(est.alpha, NA_real_, alpha),
                quietly = opt$quietly,
                input_error = prep$input_error,
                coefficient_names = "(Intercept)"
            )
        )
    }

    prev.val <- options()$warn
    options(warn = -1)
    if( est.alpha ){
        fit.args <- list(formula = yy ~ 1, phy = prep$tree, model = "BM",
                         measurement_error = opt$measurement_error)
    } else {
        fit.args <- list(formula = yy ~ 1, phy = prep$tree, model = opt$root.model,
                         measurement_error = opt$measurement_error,
                         starting.value = alpha,
                         lower.bound = alpha,
                         upper.bound = alpha)
    }
    if( !is.null(prep$input_error) ){
        fit.args$input_error <- prep$input_error
    }
    fit <- try(do.call(phylolm, fit.args), silent = opt$quietly)
    options(warn = prev.val)

    if( inherits(fit, "try-error") ){
        stop("failed to estimate observation-error parameters for whitening.")
    }

    if( is.null(fit$sigma2_error) ){
        fit$sigma2_error <- 0
    }

    return(fit)
}

extra_error_df <- function(opt){
    if( isTRUE(opt$measurement_error) ){
        return(1)
    }
    return(0)
}

phylolm_supports_input_error <- function(){
    return("input_error" %in% names(formals(phylolm)))
}

requires_phylolm_input_error_support <- function(measurement_error=FALSE, input_error=NULL){

    !is.null(input_error) &&
        isTRUE(measurement_error) &&
        !phylolm_supports_input_error()
}

can_use_dense_input_error_fit <- function(measurement_error=FALSE, input_error=NULL){

    !is.null(input_error) &&
        isFALSE(measurement_error) &&
        !phylolm_supports_input_error()
}

check_input_error_support <- function(measurement_error=FALSE, input_error=NULL){

    if( requires_phylolm_input_error_support(measurement_error, input_error) ){
        stop("input_error is not supported when measurement_error=TRUE in the current kfl1ou likelihood solver. Set measurement_error=FALSE.")
    }

    invisible(NULL)
}

dense_known_input_error_gls_fit <- function(tree, Y, preds, model,
                                            lower.bound=NA_real_,
                                            upper.bound=NA_real_,
                                            starting.value=NA_real_,
                                            quietly=TRUE,
                                            input_error,
                                            coefficient_names=NULL){

    Y <- as.matrix(Y)
    preds <- as.matrix(preds)
    y <- as.numeric(Y[, 1])
    n <- length(y)
    d <- ncol(preds)

    if( n <= d ){
        stop("not enough taxa remain to fit the model with the requested input_error.")
    }

    phy <- reorder(tree, "pruningwise")
    mean.tip.height <- mean(pruningwise.distFromRoot(phy)[seq_len(n)])

    fit_for_alpha <- function(alpha){

        root.model <- ifelse(model == "BM", "OUfixedRoot", model)
        re <- sqrt_OU_covariance(
            tree,
            alpha = ifelse(model == "BM", 0, alpha),
            root.model = root.model,
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
            if( is.null(invXX) ){
                return(NULL)
            }
            list(invXX = invXX, betahat = drop(invXX %*% Xy))
        })

        if( is.null(inv.solve) ){
            return(list(n2llh = .Machine$double.xmax))
        }

        fitted.values <- drop(preds %*% inv.solve$betahat)
        residuals <- y - fitted.values
        sigma2hat <- sum((yt - drop(Xt %*% inv.solve$betahat))^2) / n
        if( !is.finite(sigma2hat) || sigma2hat <= 0 ){
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
            residuals = residuals
        )
    }

    tol <- 1e-10
    if( model == "BM" ){
        alpha.hat <- NULL
        fit <- fit_for_alpha(0)
        p <- d + 1L
    } else{
        lower.bound <- ifelse(is.null(lower.bound), NA_real_, as.numeric(lower.bound[[1]]))
        upper.bound <- ifelse(is.null(upper.bound), NA_real_, as.numeric(upper.bound[[1]]))
        starting.value <- ifelse(is.null(starting.value), NA_real_, as.numeric(starting.value[[1]]))

        if( is.na(lower.bound) && is.na(starting.value) ){
            lower.opt <- 1e-07 / mean.tip.height
        } else{
            lower.opt <- ifelse(is.na(lower.bound), 0, lower.bound)
        }

        if( is.na(upper.bound) || upper.bound <= 0 ){
            stop("input_error fallback requires a strictly positive alpha.upper bound.")
        }

        if( lower.opt == upper.bound ){
            alpha.hat <- lower.opt
        } else if( lower.opt > 0 ){
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

        if( !quietly &&
            (isTRUE(all.equal(alpha.hat, lower.opt, tolerance = tol)) ||
             isTRUE(all.equal(alpha.hat, upper.bound, tolerance = tol))) ){
            warning(paste("the estimation of alpha matches the upper/lower bound for this parameter.\n                          You may change the bounds using options \"upper.bound\" and \"lower.bound\".\n"))
        }

        fit <- fit_for_alpha(alpha.hat)
        p <- d + 2L
    }

    if( !is.finite(fit$n2llh) || fit$n2llh >= .Machine$double.xmax ){
        stop("failed to fit the OU model with the supplied input_error.")
    }

    if( is.null(coefficient_names) ){
        coefficient_names <- paste0("preds", seq_len(d))
    }

    coefficients <- fit$betahat
    names(coefficients) <- coefficient_names
    vcov <- fit$vcov
    colnames(vcov) <- coefficient_names
    rownames(vcov) <- coefficient_names

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
        X = preds,
        n = n,
        d = d,
        model = model
    )
}

get_l1ou_thread_state <- function(){
    state <- tryCatch(get_runtime_thread_settings(),
                      error = function(e) list(openblas_threads=NA_integer_,
                                               openmp_threads=NA_integer_))
    env.names <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                   "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                   "BLIS_NUM_THREADS")
    state$env <- Sys.getenv(env.names, unset = NA_character_)
    return(state)
}

set_l1ou_thread_limit <- function(nThreads){

    nThreads <- max(1L, as.integer(nThreads[[1]]))
    state <- get_l1ou_thread_state()

    env.vars <- setNames(rep(as.character(nThreads), length(state$env)),
                         names(state$env))
    do.call(Sys.setenv, as.list(env.vars))
    try(set_runtime_thread_settings(nThreads, nThreads), silent=TRUE)

    return(state)
}

restore_l1ou_thread_limit <- function(state){

    env.state <- state$env
    if(length(env.state) > 0){
        unset <- names(env.state)[is.na(env.state)]
        if(length(unset) > 0){
            Sys.unsetenv(unset)
        }
        restore <- env.state[!is.na(env.state)]
        if(length(restore) > 0){
            do.call(Sys.setenv, as.list(restore))
        }
    }

    blas.threads <- ifelse(is.null(state$openblas_threads), NA_integer_,
                           as.integer(state$openblas_threads))
    omp.threads <- ifelse(is.null(state$openmp_threads), NA_integer_,
                          as.integer(state$openmp_threads))
    try(set_runtime_thread_settings(
            ifelse(is.na(blas.threads), -1L, blas.threads),
            ifelse(is.na(omp.threads), -1L, omp.threads)
        ),
        silent=TRUE)

    invisible(NULL)
}

with_l1ou_thread_limit <- function(nThreads, expr){

    state <- set_l1ou_thread_limit(nThreads)
    on.exit(restore_l1ou_thread_limit(state), add=TRUE)
    force(expr)
}

resolve_l1ou_thread_limit <- function(nCores=1, parallel.computing=FALSE){

    nCores <- as.integer(nCores[[1]])
    if(is.na(nCores) || nCores < 1){
        nCores <- 1L
    }
    if(parallel.computing){
        return(1L)
    }
    return(nCores)
}

make_nested_l1ou_options <- function(opt){

    nested.opt <- opt
    nested.opt$nCores <- 1L
    nested.opt$parallel.computing <- FALSE
    nested.opt$prepared.tree <- NULL
    nested.opt$prepared.tree.list <- NULL
    return(nested.opt)
}

prepare_fast_phylolm_tree <- function(tree){

    phy <- reorder(tree, "pruningwise")
    n <- length(phy$tip.label)
    des <- phy$edge[, 2]
    dist.from.root <- pruningwise.distFromRoot(phy)
    tip.height <- mean(dist.from.root[seq_len(n)])
    times <- tip.height - dist.from.root[(n + 1):(n + phy$Nnode)]

    list(
        phy = phy,
        n = n,
        N = nrow(phy$edge),
        ROOT = n + 1L,
        anc = phy$edge[, 1],
        des = des,
        externalEdge = des <= n,
        times = times,
        edge.age = times[phy$edge[, 1] - n],
        tip.height = tip.height,
        Tmax = max(times)
    )
}

edge_scaling_from_cache <- function(edge.age, type=c("apprX", "orgX"), alpha=NA_real_){

    type <- match.arg(type)

    if(type == "orgX"){
        return(1 - exp(-alpha * edge.age))
    }
    return(edge.age)
}

cached_weighted_design_matrix <- function(Z, edge.age, type=c("apprX", "orgX"), alpha=NA_real_){

    scales <- edge_scaling_from_cache(edge.age, type=type, alpha=alpha)
    sweep(Z, 2, scales, FUN="*")
}

initialize_design_cache <- function(tree, opt){

    if(is.null(opt$Z)){
        opt$Z <- generate_design_matrix(tree, "simpX")
    }
    if(is.null(opt$edge.age)){
        opt$edge.age <- tree_edge_ages(tree)
    }
    if(!is.null(opt$tree.list)){
        opt$tree.list <- lapply(opt$tree.list, function(tr){
            if(is.null(tr$Z)){
                tr$Z <- generate_design_matrix(tr, "simpX")
            }
            if(is.null(tr$edge.age)){
                tr$edge.age <- tree_edge_ages(tr)
            }
            tr
        })
    }

    return(opt)
}

initialize_fast_phylolm_cache <- function(tree, opt){

    opt$prepared.tree <- prepare_fast_phylolm_tree(tree)
    if(!is.null(opt$tree.list)){
        opt$prepared.tree.list <- lapply(opt$tree.list, prepare_fast_phylolm_tree)
    } else{
        opt$prepared.tree.list <- NULL
    }
    return(opt)
}

fast_phylolm_ou_fit <- function(prepared.tree, Y, preds, opt){

    phy <- prepared.tree$phy
    n <- prepared.tree$n
    N <- prepared.tree$N
    d <- ncol(preds)
    y <- as.numeric(Y[, 1])
    X <- as.matrix(preds)
    ole <- 4 + 2*d + d*d
    tol <- 1e-10
    loglik.cache <- new.env(parent=emptyenv())
    loglik.cache$key <- NULL
    loglik.cache$value <- NULL

    loglik <- function(alpha){
        cache.key <- format(alpha, digits=17, scientific=FALSE, trim=TRUE)
        if( !is.null(loglik.cache$key) && identical(loglik.cache$key, cache.key) ){
            return(loglik.cache$value)
        }

        if(opt$root.model == "OUrandomRoot"){
            distFromRoot <- exp(-2 * alpha * prepared.tree$times)
            d1 <- distFromRoot[prepared.tree$anc - n]
            d2 <- numeric(N)
            d2[prepared.tree$externalEdge] <- 1
            d2[!prepared.tree$externalEdge] <- distFromRoot[
                prepared.tree$des[!prepared.tree$externalEdge] - n
            ]
        } else{
            distFromRoot <- exp(-2 * alpha * prepared.tree$times) *
                (1 - exp(-2 * alpha * (prepared.tree$Tmax - prepared.tree$times)))
            d1 <- distFromRoot[prepared.tree$anc - n]
            d2 <- numeric(N)
            d2[prepared.tree$externalEdge] <- 1 - exp(-2 * alpha * prepared.tree$Tmax)
            d2[!prepared.tree$externalEdge] <- distFromRoot[
                prepared.tree$des[!prepared.tree$externalEdge] - n
            ]
        }

        edge.length <- d2 - d1
        tmp <- threepoint_l1ou_c(
            as.integer(N),
            as.integer(n),
            as.integer(phy$Nnode),
            as.integer(1),
            as.integer(d),
            as.integer(prepared.tree$ROOT),
            as.double(min(distFromRoot)),
            as.double(edge.length),
            as.integer(prepared.tree$des),
            as.integer(prepared.tree$anc),
            as.double(y),
            as.double(as.vector(X))
        )

        XX <- matrix(tmp[(5 + d):(ole - d)], d, d)
        Xy <- tmp[(ole - d + 1):ole]
        inv.solve <- tryCatch({
            chol.XX <- chol(XX)
            list(invXX = chol2inv(chol.XX),
                 betahat = backsolve(chol.XX, forwardsolve(t(chol.XX), Xy)))
        }, error = function(e) {
            invXX <- solve(XX)
            list(invXX = invXX, betahat = drop(invXX %*% Xy))
        })
        invXX <- inv.solve$invXX
        betahat <- as.numeric(inv.solve$betahat)

        sigma2hat <- as.numeric(
            (tmp[4] - 2 * sum(betahat * Xy) + crossprod(betahat, XX %*% betahat)) / n
        )
        if(sigma2hat < 0){
            resdl <- X %*% betahat - y
            sigma2hat <- threepoint_l1ou_c(
                as.integer(N),
                as.integer(n),
                as.integer(phy$Nnode),
                as.integer(1),
                as.integer(d),
                as.integer(prepared.tree$ROOT),
                as.double(min(distFromRoot)),
                as.double(edge.length),
                as.integer(prepared.tree$des),
                as.integer(prepared.tree$anc),
                as.double(as.vector(resdl)),
                as.double(as.vector(X))
            )[4] / n
        }
        vcov <- sigma2hat * invXX * n/(n - d)
        res <- list(
            n2llh = as.numeric(n * log(2 * pi) + n + n * log(sigma2hat) + tmp[1]),
            betahat = betahat,
            sigma2hat = sigma2hat,
            vcov = vcov
        )
        loglik.cache$key <- cache.key
        loglik.cache$value <- res
        res
    }

    if( is.na(opt$alpha.lower.bound) & is.na(opt$alpha.starting.value) ){
        lower <- 1e-07 / prepared.tree$tip.height
    } else{
        lower <- ifelse(is.na(opt$alpha.lower.bound), 0, opt$alpha.lower.bound)
    }
    upper <- opt$alpha.upper.bound

    if(lower == upper){
        if(lower <= 0){
            stop("fast OU fit requires a strictly positive fixed alpha.")
        }
        alpha.hat <- lower
    } else{
        if(upper <= 0){
            stop("fast OU fit requires a strictly positive upper alpha bound.")
        }
        lower.opt <- ifelse(lower <= 0, 1e-07 / prepared.tree$tip.height, lower)
        if(lower.opt >= upper){
            alpha.hat <- lower.opt
        } else{
            ## The fast path only optimizes a single positive parameter, so
            ## Brent's method is substantially cheaper than a generic L-BFGS-B run.
            opt.res <- optimize(
                function(logalpha) loglik(exp(logalpha))$n2llh,
                interval = c(log(lower.opt), log(upper)),
                tol = 1e-6
            )
            alpha.hat <- as.numeric(exp(opt.res$minimum))
        }
    }

    if( !opt$quietly &&
        (isTRUE(all.equal(alpha.hat, lower, tol = tol)) ||
         isTRUE(all.equal(alpha.hat, upper, tol = tol))) ){
        warning(paste("the estimation of alpha matches the upper/lower bound for this parameter.\n                          You may change the bounds using options \"upper.bound\" and \"lower.bound\".\n"))
    }

    fit <- loglik(alpha.hat)
    sigma2hat <- 2 * alpha.hat * fit$sigma2hat
    coefficients <- fit$betahat
    names(coefficients) <- colnames(X)
    vcov <- fit$vcov
    colnames(vcov) <- colnames(X)
    rownames(vcov) <- colnames(X)

    results <- list(
        coefficients = coefficients,
        sigma2 = sigma2hat,
        optpar = alpha.hat,
        sigma2_error = 0,
        logLik = -fit$n2llh/2,
        p = 2 + d,
        aic = 2 * (2 + d) + fit$n2llh,
        vcov = vcov,
        fitted.values = drop(X %*% coefficients),
        residuals = y - drop(X %*% coefficients),
        mean.tip.height = prepared.tree$tip.height,
        y = y,
        X = X,
        n = n,
        d = d,
        model = opt$root.model
    )
    return(results)
}

get_data <- function(tree, Y, shift.configuration, opt, idx){

    if(!is.null(opt$tree.list)){
        tr    <- opt$tree.list[[idx]]
        y.ava <- !is.na(Y[,idx])
        y     <- as.matrix(Y[y.ava, idx])
        mapped <- tr$old.order[shift.configuration]
        s.c <- as.integer(mapped[!is.na(mapped)])
        stopifnot(length(tr$tip.label)==nrow(y))
        input.error <- get_trait_input_error(opt, idx, tree=tr, available=y.ava)
        prepared.tree <- opt$prepared.tree.list[[idx]]
    }else{
        tr  <- tree
        y   <- as.matrix(Y[, idx])
        s.c <- shift.configuration
        input.error <- get_trait_input_error(opt, idx, tree=tr)
        prepared.tree <- opt$prepared.tree
    }
    result     <- list()
    result$tr  <- tr
    result$y   <- y
    result$s.c <- s.c
    result$input_error <- input.error
    result$prepared.tree <- prepared.tree
    return(result)
}

cmp_BIC <- function(tree, Y, shift.configuration, opt){

    nEdges     <- Nedge(tree)
    nTips      <- length(tree$tip.label)
    nShifts    <- length(shift.configuration)
    nVariables <- ncol(Y)

    df.1  <- log(nTips)*(nShifts)
    score <- df.1
    alpha <- sigma2 <- logLik <- rep(0, nVariables)

    for( i in 1:nVariables ){

        r   <- get_data(tree, Y, shift.configuration, opt, i)
        tr  <- r$tr
        y   <- r$y
        s.c <- r$s.c

        fit  <- my_phylolm_interface(tr, y, s.c, opt, input_error=r$input_error,
                                     prepared.tree = r$prepared.tree)
        if ( all(is.na(fit)) ){ return(NA) } 

        df.2 <- log(nrow(y))*fit$p
        score <- score  -2*fit$logLik + df.2

        alpha [[i]] <- fit$optpar
        sigma2[[i]] <- fit$sigma2
        logLik[[i]] <- fit$logLik
    }
    return( list(score=score, alpha=alpha, sigma2=sigma2, logLik=logLik) )
}


cmp_AICc <- function(tree, Y, shift.configuration, opt){

    nShifts    <- length(shift.configuration)
    nVariables <- ncol(Y)

    score <- 0
    total.p <- 0
    total.n <- 0
    alpha <- sigma2 <- logLik <- rep(0, nVariables)

    for( i in 1:nVariables ){

        r   <- get_data(tree, Y, shift.configuration, opt, i)
        tr  <- r$tr
        y   <- r$y
        s.c <- r$s.c

        fit <- my_phylolm_interface(tr, y, s.c, opt, input_error=r$input_error,
                                    prepared.tree = r$prepared.tree)
        if ( all(is.na(fit)) ){ return(NA) } 
        score <- score  -2*fit$logLik
        total.p <- total.p + fit$p
        total.n <- total.n + fit$n

        alpha [[i]] <- fit$optpar
        sigma2[[i]] <- fit$sigma2
        logLik[[i]] <- fit$logLik
    }

    p <- nShifts + total.p
    if( p > total.n - 2 )
        return(NA)

    d.f <- 2*p + (2*p*(p+1))/(total.n-p-1)
    score <- score + d.f

    return( list(score=score, alpha=alpha, sigma2=sigma2, logLik=logLik) )
}

cmp_mBIC <- function(tree, Y, shift.configuration, opt){

    nEdges     <- Nedge(tree)
    nTips      <- length(tree$tip.label)
    nShifts    <- length(shift.configuration)
    nVariables <- ncol(Y)

    res =  cmp_mBIC_df(tree, shift.configuration, opt)  
    df.1 = res$df.1
    df.2 = res$df.2

    score <- df.1

    alpha <- sigma2 <- logLik <- rep(0, nVariables)
    for( i in 1:nVariables ){

        r   <- get_data(tree, Y, shift.configuration, opt, i)
        tr  <- r$tr
        y   <- r$y
        s.c <- r$s.c

        if( nVariables > 1){
            res  = cmp_mBIC_df(tr, s.c, opt)  
            df.2 = res$df.2
        }

        fit <- my_phylolm_interface(tr, y, s.c, opt, input_error=r$input_error,
                                    prepared.tree = r$prepared.tree)
        if ( all(is.na(fit)) ){ return(NA) } 

        score <- score  -2*fit$logLik + df.2

        alpha [[i]] <- fit$optpar
        sigma2[[i]] <- fit$sigma2
        logLik[[i]] <- fit$logLik
    }
    return( list(score=score, alpha=alpha, sigma2=sigma2, logLik=logLik) )
}

cmp_mBIC_df <- function(tree, shift.configuration, opt){

    if(length(shift.configuration)>0){
        shift.configuration <- sort(shift.configuration)
    }
    nTips               <- length(tree$tip.label)
    nShifts             <- length(shift.configuration)

    df.1 <- 0 
    ## pen for the alpha sigma2 and intercept
    df.2 <- (3 + extra_error_df(opt))*log(nTips)

    if(nShifts > 0 ){
        ## pen for shift configuration
        df.1 = (2*nShifts - 1) *log(nTips)
        ## pen for alpha sigma2 and intercept
        df.2 = (3 + extra_error_df(opt))*log(nTips)

        all.covered.tips = rep(FALSE, nTips)
        for(eIdx in shift.configuration){
            covered.tips = opt$Z[,eIdx] > 0
            nUniqueTips  = sum(covered.tips & !all.covered.tips)
            all.covered.tips = all.covered.tips | covered.tips
            ## this must not happen if the input is an 
            ## identifiable configuration (parsimonious) and the tree is in post order.
            stopifnot( nUniqueTips > 0)
            df.2 = df.2 + log(nUniqueTips) 
        }
        nUniqueTips = sum(!all.covered.tips)
        df.2 = df.2 + log(nUniqueTips) 
    } 

    return( list(df.1=df.1, df.2=df.2) )
}

cmp_pBICess <- function(tree, Y, shift.configuration, opt){

    nShifts = length(shift.configuration)
    nEdges  = Nedge(tree)
    nTips   = length(tree$tip.label)

    df.1  = 2*(nShifts)*log(nEdges-1)
    score = df.1

    alpha = sigma2  = logLik = numeric()

    for(i in 1:ncol(Y)){

        r   <- get_data(tree, Y, shift.configuration, opt, i)
        tr  <- r$tr
        y   <- r$y
        s.c <- r$s.c

        fit  = my_phylolm_interface(tr, y, s.c, opt, input_error=r$input_error,
                                    prepared.tree = r$prepared.tree)
        if( all(is.na(fit)) ){
           return(NA)
        }
        ess  = effective.sample.size(tr, edges=s.c, model="OUfixedRoot", 
                 parameters=list(alpha=fit$optpar), FALSE, FALSE)

        df.2  = (3 + extra_error_df(opt))*log(nrow(y)+1) + sum(log(ess+1))
        score = score  -2*fit$logLik + df.2 

        alpha  = c(alpha, fit$optpar)
        sigma2 = c(sigma2, fit$sigma2)
        logLik = c(logLik, fit$logLik)
    }
    return( list(score=score, alpha=alpha, sigma2=sigma2, logLik=logLik) )
}


cmp_pBIC <- function(tree, Y, shift.configuration, opt){

    nShifts = length(shift.configuration)
    nEdges  = Nedge(tree)
    nTips   = length(tree$tip.label)

    df.1    = 2*(nShifts)*log(nEdges-1)
    score   = df.1
    alpha   = sigma2  = logLik = rep(0, ncol(Y))

    for(i in 1:ncol(Y)){

        r   <- get_data(tree, Y, shift.configuration, opt, i)
        tr  <- r$tr
        y   <- r$y
        s.c <- r$s.c

        fit   = my_phylolm_interface(tr, y, s.c, opt, input_error=r$input_error,
                                     prepared.tree = r$prepared.tree)
        if( all(is.na(fit)) ){
           return(NA)
        } 
        varY  = c(var(y))
        ld    = as.numeric(determinant(fit$vcov * (fit$n - fit$d)/(varY*fit$n),
                                      logarithm = TRUE)$modulus)
        df.2  = (2 + extra_error_df(opt))*log(nrow(y)) - ld
        score = score  -2*fit$logLik + df.2 

        alpha [[i]] = fit$optpar
        sigma2[[i]] = fit$sigma2
        logLik[[i]] = fit$logLik
    }
    return( list(score=score, alpha=alpha, sigma2=sigma2, logLik=logLik) )
}

my_phylolm_interface <- function(tree, Y, shift.configuration, opt, recmp.preds=FALSE,
                                 alpha=0, input_error=NULL, prepared.tree=NULL){

    if(recmp.preds){
        Z <- generate_design_matrix(tree, type="orgX", alpha=alpha)
    }else{
        if(opt$multivariate.missing){
            if(is.null(tree$Z))
                stop("internal error: in multivariate.missing mode but tree$Z is null.")
            Z <- tree$Z
        }else{
            Z <- opt$Z 
        }
    }


    preds = cbind(1, Z[ ,shift.configuration])
    if( requires_phylolm_input_error_support(opt$measurement_error, input_error) ){
        stop("input_error is not supported when measurement_error=TRUE in the current kfl1ou likelihood solver. Set measurement_error=FALSE.")
    }
    if( can_use_dense_input_error_fit(opt$measurement_error, input_error) ){
        fit <- try(
            dense_known_input_error_gls_fit(
                tree,
                Y,
                preds,
                model = opt$root.model,
                lower.bound = opt$alpha.lower.bound,
                upper.bound = opt$alpha.upper.bound,
                starting.value = opt$alpha.starting.value,
                quietly = opt$quietly,
                input_error = input_error,
                coefficient_names = paste0("preds", seq_len(ncol(preds)))
            ),
            silent = opt$quietly
        )
        if(class(fit) != "try-error"){
            return(fit)
        }
    }

    if( isFALSE(opt$measurement_error) &&
        is.null(input_error) &&
        !recmp.preds &&
        opt$root.model %in% c("OUfixedRoot", "OUrandomRoot") &&
        !is.null(prepared.tree) ){
        fit <- try(fast_phylolm_ou_fit(prepared.tree, Y, preds, opt),
                   silent = opt$quietly)
        if(class(fit) != "try-error"){
            return(fit)
        }
    }

    prev.val <-options()$warn 
    options(warn = -1)
    if( is.na(opt$alpha.lower.bound) & is.na(opt$alpha.starting.value) ){
        fit.args <- list(formula = Y~preds-1, phy=tree, model=opt$root.model,
                         upper.bound  = opt$alpha.upper.bound,
                         measurement_error = opt$measurement_error)
    } else {
        l = ifelse(is.na(opt$alpha.lower.bound), 0, opt$alpha.lower.bound)
        u = opt$alpha.upper.bound
        s = ifelse(is.na(opt$alpha.starting.value), max(0.5, l), opt$alpha.starting.value)

        fit.args <- list(formula = Y~preds-1, phy = tree,
                         model          = opt$root.model,
                         starting.value = s,
                         lower.bound    = l,
                         upper.bound    = u,
                         measurement_error = opt$measurement_error)
    }
    if( !is.null(input_error) ){
        fit.args$input_error <- input_error
    }
    fit <- try(do.call(phylolm, fit.args), silent = opt$quietly)
    options(warn = prev.val )

    if(class(fit) == "try-error"){ 
        if(!opt$quietly){
            warning( paste0( "the OU likelihood solver returned an error with a shift configuration
                            of size ", length(shift.configuration), ". You may
                            want to change alpha.upper/alpha.lower!") )
        }
        return(NA)
    }
    return(fit)
}


run_grplasso  <- function (grpX, grpY, nVariables, grpIdx, opt){
    delta  = opt$grp.delta
    seq.ub = opt$grp.seq.ub
    backend <- resolve_grplasso_backend(opt)
    max.nTries = 7
    coarse.delta <- max(delta * 4, delta)
    min.delta <- delta / 4
    lmbdMax  <-  1.2 * linreg_group_lasso_lambda_max(grpX, grpY, grpIdx, backend = backend) + 1

    base.seq = seq(0, seq.ub, coarse.delta)
    for (itrTmp in 1:max.nTries) {
        lmbd = lmbdMax * (0.5^base.seq)
        sol <- run_grplasso_path(grpX, grpY, grpIdx, lmbd, tol = 0.01, backend = backend)

        df.vec <- count_active_grplasso_groups(sol$coefficients, grpIdx, nVariables)

        next.seq <- grplasso_refine_base_seq(base.seq, df.vec, opt$max.nShifts, min.delta)
        if(length(next.seq) == length(base.seq) &&
           isTRUE(all.equal(next.seq, base.seq, tolerance = 1e-12))){
            indices = which(df.vec > (opt$max.nShifts + 4))
            if (length(indices) > 0) {
                upper.idx = min(indices)
                base.seq  = base.seq[1:upper.idx]
            }
            break 
        } else{
            base.seq = next.seq
        }
    }

    final.lmbd <- lmbdMax * (0.5^base.seq)
    reuse.coarse <- exists("sol") && exists("lmbd") &&
        length(lmbd) >= length(final.lmbd) &&
        isTRUE(all.equal(as.numeric(lmbd[seq_along(final.lmbd)]),
                         as.numeric(final.lmbd), tolerance = 1e-12))
    if(reuse.coarse){
        coarse.sol <- sol
        coarse.sol$coefficients <- as.matrix(sol$coefficients)[, seq_along(final.lmbd), drop=FALSE]
    } else{
        coarse.sol <- run_grplasso_path(grpX, grpY, grpIdx, final.lmbd, tol = 0.01, backend = backend)
    }
    lmbd <- final.lmbd
    df.vec <- count_active_grplasso_groups(coarse.sol$coefficients, grpIdx, nVariables)
    df.missing = setdiff(0:opt$max.nShifts, df.vec)
    refine.idx <- grplasso_support_change_indices(coarse.sol$coefficients, grpIdx, nVariables)
    if(length(refine.idx) == 0){
        refine.idx <- 1L
    }
    sol <- run_grplasso_path(grpX, grpY, grpIdx, lmbd[refine.idx], tol = 1e-6, backend = backend)

    for (dfm in df.missing) {
        warning(paste0("There are no solutions with ", dfm, " number of shifts  
                in the solution path of grplasso. You may want to change grp.delta and grp.seq"))
    }
    return(sol);
}

grplasso_group_nonzero_counts <- function(coefficients, grpIdx){

    valid <- !is.na(grpIdx)
    if(!any(valid)){
        return(matrix(0L, nrow=0, ncol=ncol(as.matrix(coefficients))))
    }

    coeff.mat <- as.matrix(coefficients)[valid, , drop=FALSE]
    rowsum((abs(coeff.mat) > 0) + 0L, group=grpIdx[valid], reorder=FALSE)
}

count_active_grplasso_groups <- function(coefficients, grpIdx, nVariables){

    active.by.group <- grplasso_group_nonzero_counts(coefficients, grpIdx)
    if(nrow(active.by.group) == 0){
        return(integer(ncol(as.matrix(coefficients))))
    }
    threshold <- max(1L, ceiling(nVariables / 2))
    colSums(active.by.group >= threshold)
}

grplasso_support_change_indices <- function(coefficients, grpIdx, nVariables){

    nSolutions <- ncol(as.matrix(coefficients))
    if(nSolutions == 0){
        return(integer())
    }

    active.by.group <- grplasso_group_nonzero_counts(coefficients, grpIdx)
    if(nrow(active.by.group) == 0){
        return(1L)
    }

    threshold <- max(1L, ceiling(nVariables / 2))
    support.mat <- active.by.group >= threshold
    signatures <- vapply(seq_len(ncol(support.mat)), function(idx) {
        paste0(which(support.mat[, idx]), collapse=" ")
    }, character(1))
    keep <- c(TRUE, signatures[-1] != signatures[-nSolutions])
    which(keep)
}

grplasso_refine_base_seq <- function(base.seq, df.vec, max.nShifts, min.delta){

    targets <- setdiff(0:(max.nShifts + 1L), unique(df.vec))
    if(length(targets) == 0){
        return(base.seq)
    }

    if(length(base.seq) > 1){
        steps <- diff(base.seq)
    } else{
        steps <- numeric()
    }
    df.search <- cummax(df.vec)
    new.points <- numeric()

    if(length(base.seq) > 1){
        for(target in targets){
            idx <- which(
                pmin(df.search[-length(df.search)], df.search[-1]) <= target &
                pmax(df.search[-length(df.search)], df.search[-1]) >= target &
                df.search[-length(df.search)] != df.search[-1] &
                steps > min.delta
            )
            if(length(idx) > 0){
                new.points <- c(new.points, (base.seq[idx] + base.seq[idx + 1L]) / 2)
            }
        }
    }

    if(max(df.search) <= max.nShifts){
        extension.step <- if(length(steps) > 0) tail(steps, 1) else max(min.delta * 2, 1)
        if(extension.step > min.delta){
            new.points <- c(new.points, base.seq[[length(base.seq)]] + extension.step)
        }
    }

    unique(sort(c(base.seq, new.points)))
}
