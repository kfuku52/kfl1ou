#' Writes placeholder outputs for workflows expecting l1ou-style files
#'
#' This helper writes placeholder TSV and PDF outputs for workflows that expect
#' \code{l1ou_tree.tsv}, \code{l1ou_regime.tsv}, \code{l1ou_leaf.tsv}, and a
#' simple tree plot even when no shifts can be fitted, for example because all
#' traits are invariant.
#'
#'@param tree tree of class \code{phylo}.
#'@param Y trait matrix or data frame aligned to \code{tree$tip.label}.
#'@param dir output directory.
#'@param prefix output file prefix. The default \code{"l1ou"} writes files such
#'  as \code{l1ou_tree.tsv}.
#'@param plot.width width of the placeholder PDF in inches.
#'@param plot.height optional height of the placeholder PDF in inches. When
#'  \code{NULL}, the height scales with the number of tips.
#'@param reason optional character string stored in the returned metadata.
#'
#'@return Invisibly returns a list containing placeholder tables, file paths,
#'  and the optional \code{reason}.
#'
#'@examples
#' data(lizard.tree, lizard.traits)
#' keep <- lizard.tree$tip.label[1:8]
#' tree <- ape::drop.tip(lizard.tree, setdiff(lizard.tree$tip.label, keep))
#' Y <- as.matrix(lizard.traits[keep, 1, drop = FALSE])
#' Y[] <- 1
#' out.dir <- tempdir()
#' write_l1ou_placeholder_outputs(tree, Y, dir = out.dir)
#'
#'@export
write_l1ou_placeholder_outputs <- function(tree, Y, dir = ".", prefix = "l1ou",
                                           plot.width = 8, plot.height = NULL,
                                           reason = NULL){

    if( !inherits(tree, "phylo") ){
        stop("tree must be of class \"phylo\".")
    }

    trait.table <- as_placeholder_trait_table(tree, Y)

    if( is.null(plot.height) ){
        plot.height <- length(tree$tip.label) / 5 + 1
    }

    if( !dir.exists(dir) ){
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }

    tree.table <- get_l1ou_placeholder_tree_table(tree, trait.table)
    regime.table <- get_l1ou_placeholder_regime_table(tree, trait.table)
    leaf.table <- get_l1ou_placeholder_leaf_table(tree, trait.table)

    tree.file <- file.path(dir, paste0(prefix, "_tree.tsv"))
    regime.file <- file.path(dir, paste0(prefix, "_regime.tsv"))
    leaf.file <- file.path(dir, paste0(prefix, "_leaf.tsv"))
    plot.file <- file.path(dir, paste0(prefix, "_plot.pdf"))

    write.table(tree.table, file = tree.file, sep = "\t", quote = FALSE,
                row.names = FALSE)
    write.table(regime.table, file = regime.file, sep = "\t", quote = FALSE,
                row.names = FALSE)
    write.table(leaf.table, file = leaf.file, sep = "\t", quote = FALSE,
                row.names = FALSE)

    grDevices::pdf(plot.file, height = plot.height, width = plot.width)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot(tree, show.tip.label = TRUE, cex = 0.5)

    invisible(list(
        tree_table = tree.table,
        regime_table = regime.table,
        leaf_table = leaf.table,
        files = c(tree = tree.file,
                  regime = regime.file,
                  leaf = leaf.file,
                  plot = plot.file),
        reason = reason
    ))
}

as_placeholder_trait_table <- function(tree, Y){

    if( inherits(Y, "data.frame") ){
        Y <- as.data.frame(Y, check.names = FALSE, stringsAsFactors = FALSE)
    } else{
        Y <- as.data.frame(as.matrix(Y), check.names = FALSE, stringsAsFactors = FALSE)
    }

    if( nrow(Y) != length(tree$tip.label) ){
        stop("Y must have one row per tip in tree.")
    }

    if( is.null(rownames(Y)) ){
        rownames(Y) <- tree$tip.label
    }

    if( !setequal(rownames(Y), tree$tip.label) ){
        stop("rownames(Y) must match tree$tip.label.")
    }

    Y <- Y[tree$tip.label, , drop = FALSE]
    return(Y)
}

get_l1ou_placeholder_leaf_table <- function(tree, original_trait_table){

    out <- data.frame()
    params <- c("Y", "optima", "mu", "residuals")
    for(param in params){
        tmp <- data.frame(
            regime = rep(0, nrow(original_trait_table)),
            node_name = rownames(original_trait_table),
            param = param,
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
        tmp <- cbind(tmp, original_trait_table)
        rownames(tmp) <- NULL
        out <- rbind(out, tmp)
    }
    return(out)
}

get_l1ou_placeholder_regime_table <- function(tree, original_trait_table){

    out <- data.frame()
    params <- c("alpha", "sigma2", "intercept", "log_likelihood")
    trait.cols <- colnames(original_trait_table)
    if( is.null(trait.cols) ){
        trait.cols <- paste0("trait", seq_len(ncol(original_trait_table)))
    }
    for(param in params){
        tmp <- c(NA, NA, param, rep(NA, ncol(original_trait_table)))
        out <- rbind(out, tmp)
    }
    colnames(out) <- c("regime", "node_name", "param", trait.cols)
    return(out)
}

get_l1ou_placeholder_tree_table <- function(tree, original_trait_table){

    out <- data.frame()
    params <- c("num_shift", "num_regime", "num_conv_regime",
                "num_uniq_regime", "num_species", "num_leaf", "model_score")
    out <- rbind(out, rep(0, length(params)))
    colnames(out) <- params
    return(out)
}
