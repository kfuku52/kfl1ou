#
#' (inverse) square root of the phylogenetic covariance
#'
#' Computes an inverse square root and square root of the phylogenetic covariance matrix,
#' under the Brownian motion (BM) or the Ornstein-Uhlenbeck (OU) model.
#' The algorithm traverses the tree only once, hence the algorithm is very fast
#' and can be applied to very big trees.
#'
#'@param tree tree of class phylo with branch lengths. If alpha>0, i.e. under the OU model, the tree has to be ultrametric.
#'@param alpha adaptation rate for the OU model. The default is 0, which corresponds to the BM mode with a fixed ancestral state at the root.
#'@param root.model ancestral state model at the root.
#'@param sigma2 phylogenetic variance rate used to standardize observation-level error terms.
#'@param sigma2_error scalar observation-level error variance to add to each tip.
#'@param input_error optional vector of tip-specific observation error variances.
#'@param check.order retained for backward compatibility. Trees are always
#' reordered to postorder before the native traversal.
#'@param check.ultrametric retained for backward compatibility. Ultrametricity
#' is always enforced when \code{alpha > 0} because otherwise the OU
#' transformation is invalid.
#'
#'@details The stationary-root covariance becomes intrinsically
#' ill-conditioned as \code{alpha} approaches zero. A warning is emitted when
#' \code{alpha} is too small relative to the tree height for reliable
#' double-precision inversion.
#'
#'@return
#' \item{sqrtInvSigma}{inverse square root of the phylogenetic covariance matrix.}
#' \item{sqrtSigma}{square root of the phylogenetic covariance matrix.}
#'
#'@examples
#'
#' data(lizard.tree)
#' example_tree <- drop.tip(
#'   lizard.tree,
#'   lizard.tree$tip.label[-seq_len(15)]
#' )
#' res <- sqrt_OU_covariance(example_tree) # alpha not provided: so BM model.
#' Sigma <- vcv(example_tree)
#' dimnames(Sigma) <- NULL
#' all.equal(res$sqrtSigma %*% t(res$sqrtSigma), Sigma) # TRUE
#' all.equal(res$sqrtInvSigma %*% t(res$sqrtInvSigma), solve(Sigma)) # TRUE
#' 
#' 
#' ##Here's the example from "Eric A. Stone. 2011." (See references)
#'
#' tr <-  read.tree(text="((((Homo:.21,Pongo:.21):.28,Macaca:.49):.13,Ateles:.62):.38,Galago:1);") 
#' RE <- sqrt_OU_covariance(tr) 
#' B <- round( RE$sqrtSigma, digits=3)
#' D <- round( RE$sqrtInvSigma, digits=3)
#' print(B)
#' print(D)
#' 
#' 
#' ##Here is the examples on how to get the contrasts using sqrt_OU_covariance
#' data(lizard.tree, lizard.traits)
#' keep <- lizard.tree$tip.label[1:15]
#' tree <- drop.tip(lizard.tree, setdiff(lizard.tree$tip.label, keep))
#' tree <- reorder(tree, "postorder")
#' lizard <- adjust_data(tree, lizard.traits[keep, 1])
#' eModel <- fit_OU(lizard$tree, lizard$Y, shift.configuration=1L, criterion="AICc")
#' theta <- eModel$intercept + convert_shifts2regions(eModel$tree,
#'                              eModel$shift.configuration, eModel$shift.values)
#' REf <- sqrt_OU_covariance(eModel$tree, alpha=eModel$alpha,
#'                                          root.model = "OUfixedRoot",
#'                                          check.order=FALSE, check.ultrametric=FALSE)
#'  covInverseSqrtf  <- t(REf$sqrtInvSigma)
#'  covSqrtf   <- REf$sqrtSigma
#' # `covInverseSqrtf` is the transposed square root of the inverse
#' # covariance matrix for the fixed-root model.
#' # `covSqrtf` represents the square root of the covariance matrix for FixedRoot model.
#'  Y  <- rTraitCont(eModel$tree, "OU", theta=theta,
#'                                      alpha=eModel$alpha,
#'                                      sigma=sqrt(eModel$sigma2), root.value=eModel$intercept)
#'  contrast    <-  covInverseSqrtf%*%(Y - eModel$mu)
#'  
#'
#'@references
#' Mohammad Khabbazian, Ricardo Kriebel, Karl Rohe, and Cécile Ané (2016).
#' "Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models".
#' Methods in Ecology and Evolution. doi:10.1111/2041-210X.12534
#'
#' Eric A. Stone. 2011. "Why the phylogenetic regression appears robust to tree misspecification". Systematic Biology, 60(3):245-260.
#'
#'@export
sqrt_OU_covariance <- function(tree, alpha=0, root.model = c("OUfixedRoot", "OUrandomRoot"),
                               sigma2 = 1, sigma2_error = 0, input_error = NULL,
                               check.order=TRUE, check.ultrametric=TRUE){
    validate_l1ou_tree(tree, require.positive.edges=FALSE)
    if(length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) || alpha < 0){
        stop("alpha must be a single finite non-negative number.")
    }
    if(length(sigma2) != 1L || !is.numeric(sigma2) || !is.finite(sigma2) ||
       sigma2 <= 0){
        stop("sigma2 must be strictly positive.")
    }
    if(length(sigma2_error) != 1L || !is.numeric(sigma2_error) ||
       !is.finite(sigma2_error)){
        stop("sigma2_error must be a single finite number.")
    }
    child.count <- table(tree$edge[, 1L])
    if(any(child.count != 2L)){
        tree         <- multi2di(tree, random=FALSE)
        check.order  <- TRUE 
    }
    root.model <- match.arg(root.model) 

    if( sigma2_error < 0 ){
        stop("sigma2_error must be non-negative.")
    }
    if( !is.null(input_error) ){
        if( length(input_error) != length(tree$tip.label) ){
            stop("the number of input errors does not match the number of tips in the tree.")
        }
        if( is.null(names(input_error)) ){
            names(input_error) <- tree$tip.label
        } else {
            if(anyNA(names(input_error)) || any(!nzchar(names(input_error))) ||
               anyDuplicated(names(input_error))){
                stop("input_error names must be non-missing, non-empty, and unique.")
            }
            o <- match(tree$tip.label, names(input_error))
            if( any(is.na(o)) ){
                stop("input_error names do not match the tree tip labels.")
            }
            input_error <- input_error[o]
        }
        if(!is.numeric(input_error)){
            stop("input_error must be numeric.")
        }
        if(anyNA(input_error)){
            stop("input_error cannot contain missing values.")
        }
        if(any(!is.finite(input_error)) || any(input_error < 0)){
            stop("input_error must contain finite non-negative values.")
        }
    }
    if( (sigma2_error > 0 || !is.null(input_error)) && sigma2 <= 0 ){
        stop("sigma2 must be strictly positive when observation error is provided.")
    }

    ## Reordering changes edge rows but preserves node numbers, which the
    ## native traversal uses to associate columns with internal nodes.
    tree <- reorder(tree, "postorder")

    if ( alpha > 0){
        ## This transformation is valid only for ultrametric trees.
        if(!isTRUE(is.ultrametric(tree))){
            stop("alpha>0, the tree has to be ultrametric")
        }
        if(root.model == "OUrandomRoot"){
            tip.height <- max(
                tree_node_depths(tree)[seq_along(tree$tip.label)]
            )
            if(alpha * tip.height <= sqrt(.Machine$double.eps)){
                warning(
                    "OUrandomRoot covariance is numerically ill-conditioned ",
                    "because alpha is extremely small relative to tree height."
                )
            }
        }
        tre <- transf.branch.lengths(tree, model=root.model, parameters=list(alpha=alpha), check.pruningwise=F)$tree
	coe = 2*alpha 
        tre$edge.length=tre$edge.length/coe
	if(!is.null(tre$root.edge)) tre$root.edge=tre$root.edge/coe
    }else{
        tre <- tree
        if( root.model == "OUrandomRoot"){
            warning("alpha=0, BM model, the ancestral state model is changed to the OUfixedRoot")
        }
    }

    if( sigma2_error > 0 || !is.null(input_error) ){
        nTips <- length(tre$tip.label)
        externalEdge <- tre$edge[, 2] <= nTips
        total.error <- rep(sigma2_error/sigma2, nTips)
        if( !is.null(input_error) ){
            total.error <- total.error + input_error/sigma2
        }
        tre$edge.length[externalEdge] <-
            tre$edge.length[externalEdge] + total.error[tre$edge[externalEdge, 2]]
    }

    my.edge.list <- cbind(tre$edge-1, tre$edge.length) 
    tre$root.edge <- ifelse(is.null(tre$root.edge), 0, tre$root.edge)
    result       <- cmp_sqrt_OU_covariance(my.edge.list, length(tre$tip.label), tre$root.edge)
    return(result)
}
