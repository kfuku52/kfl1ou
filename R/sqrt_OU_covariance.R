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
#'@param check.order logical. If TRUE, the order will be checked to be in postorder traversal.
#'@param check.ultrametric logical. If TRUE, the tree will be checked to ultrametric.
#'
#'@return 
#' \item{sqrtInvSigma}{inverse square root of the phylogenetic covariance matrix.}
#' \item{sqrtSigma}{square root of the phylogenetic covariance matrix.}
#'
#'@examples
#'
#' data(lizard.tree)
#' res <- sqrt_OU_covariance(lizard.tree) # alpha not provided: so BM model.
#' Sigma <- vcv(lizard.tree)
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
#' # `covInverseSqrtf` represents the transpose of square root of  the inverse matrix of covariance for FixedRoot model.
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
    if( ! is.binary(tree) ){
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
            o <- match(tree$tip.label, names(input_error))
            if( any(is.na(o)) ){
                stop("input_error names do not match the tree tip labels.")
            }
            input_error <- input_error[o]
        }
        if( any(is.na(input_error)) ){
            stop("input_error cannot contain missing values.")
        }
        if( any(input_error < 0) ){
            stop("input_error must be non-negative.")
        }
    }
    if( (sigma2_error > 0 || !is.null(input_error)) && sigma2 <= 0 ){
        stop("sigma2 must be strictly positive when observation error is provided.")
    }

    ##NOTE: the function assumes reordering does not change the order of the 
    ##nodes and it just change the order of edges, so that column i in each 
    ##matrix still corresponds to internal node 
    ##NOTE:  in case the tree is not binary; the order will change and it is no longer postorder. 
    if( check.order ){
        tree <- reorder(tree, "post")
    }

    if ( alpha > 0){
        ##NOTE: this step requires that the tree be ultrametric tree. 
        ##NOTE: If the tree is not ultrametric, the function returns a wrong result with no warning
        if(check.ultrametric){
            if(!is.ultrametric(tree)){
                stop("alpha>0, the tree has to be ultrametric") 
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
