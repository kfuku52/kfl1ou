# Run with an installed kfl1ou: Rscript tools/reproduce-pbic-coordinates.R [csv]
# R_LIBS can select an isolated installation when comparing before and after.
suppressPackageStartupMessages(library(kfl1ou))
options(digits=17)
tree <- ape::reorder.phylo(ape::read.tree(text=paste0(
    "(((t0:1,t1:1):1,(t2:1,t3:1):1):1,",
    "((t4:1,t5:1):1,(t6:1,t7:1):1):1);"
)), "postorder")
values <- c(1.5549029785536523, 1.3137883823887941,
            -.21225628554971665, -.24002282307176648,
            1.7035058708492061, 1.8521601652123756,
            .7166573577624096, .2547434554367455)
Y <- matrix(values, ncol=1L, dimnames=list(paste0("t", 0:7), "value"))
Y <- Y[tree$tip.label, , drop=FALSE]
shifts <- match(match(c("t0", "t4"), tree$tip.label), tree$edge[, 2L])

compare_pair <- function(fixed){
    fit <- function(regimes=NULL){
        suppressWarnings(fit_OU(
            tree, Y, shifts, cr.regimes=regimes, criterion="pBIC",
            alpha.lower=if(fixed) .4 else 1e-7/3,
            alpha.upper=if(fixed) .4 else log(2)/3,
            compute.hessian=FALSE, optimizer.starts=1L
        ))
    }
    free <- fit()
    singleton <- fit(as.list(c(0L, shifts)))
    data.frame(
        alpha_mode=if(fixed) "fixed" else "estimated",
        free_score=free$score,
        singleton_score=singleton$score,
        score_difference=free$score - singleton$score,
        free_log_likelihood=sum(free$logLik),
        singleton_log_likelihood=sum(singleton$logLik),
        free_alpha=unname(free$alpha),
        singleton_alpha=unname(singleton$alpha),
        # Both chosen branches have parent age 1 in this fixture.
        coordinate_correction=2 * length(shifts) * log(-expm1(-free$alpha)),
        fixed_alpha_correction=if(fixed) -log(nrow(Y)) else 0,
        backend_version=as.character(utils::packageVersion("kfl1ou"))
    )
}

result <- do.call(rbind, lapply(c(FALSE, TRUE), compare_pair))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)) utils::write.csv(result, args[[1L]], row.names=FALSE)
print(result, row.names=FALSE)
