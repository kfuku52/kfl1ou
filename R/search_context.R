merge_l1ou_profiles <- function(profiles){
    profiles <- Filter(Negate(is.null), profiles)
    configurations <- unlist(lapply(profiles, `[[`, "configurations"),
                              recursive=FALSE)
    scores <- unlist(lapply(profiles, `[[`, "scores"), use.names=FALSE)
    if(!length(scores)) return(list(scores=numeric(), configurations=list()))
    configurations <- lapply(configurations, function(x) sort(as.integer(x)))
    keys <- vapply(configurations, configuration_score_key, character(1))
    keep <- which(!duplicated(keys, fromLast=TRUE))
    keep <- keep[order(scores[keep], keys[keep])]
    list(scores=as.numeric(scores[keep]), configurations=configurations[keep])
}

record_candidate_results <- function(results, opt){
    profiles <- lapply(results, `[[`, "profile")
    if(is.environment(opt$search.state)){
        profiles <- c(list(opt$search.state$profile), profiles)
    }
    profile <- merge_l1ou_profiles(profiles)
    if(isTRUE(opt$use.saved.scores)){
        for(i in seq_along(profile$scores)){
            add_configuration_score_to_list(
                profile$configurations[[i]], profile$scores[[i]], "",
                cache=opt$score.cache
            )
        }
    }
    if(is.environment(opt$search.state)) opt$search.state$profile <- profile
    profile
}

prepare_l1ou_search <- function(tree, Y, opt){
    requested <- opt$max.nShifts
    limit <- if(is_full_trait_covariance(opt, Y) &&
                 identical(opt$covariance.regularization, "none")){
        max(0L, nrow(Y) - ncol(Y) - 1L)
    } else requested
    opt$max.nShifts <- min(requested, limit)
    opt$search.state <- new.env(parent=emptyenv())
    opt$search.state$profile <- NULL
    opt$search.state$ensemble.paths <- list()
    opt$search.state$estimability <- list(
        requested.max.nShifts=requested,
        effective.max.nShifts=opt$max.nShifts,
        reason=if(requested > limit){
            "Unregularized full covariance requires n - rank(mean design) >= p."
        } else ""
    )
    opt
}

run_l1ou_ensemble_path <- function(task, opt, traits, rows){
    error <- ""
    result <- tryCatch(suppressWarnings(task()), error=function(e){
        error <<- conditionMessage(e)
        NULL
    })
    if(is.null(result) && !nzchar(error)) error <- "sparse path returned no solution"
    if(is.environment(opt$search.state)){
        paths <- opt$search.state$ensemble.paths
        paths[[length(paths) + 1L]] <- list(
            traits=traits, observations=length(rows),
            successful=!is.null(result), error=error
        )
        opt$search.state$ensemble.paths <- paths
    }
    result
}

l1ou_search_details <- function(opt){
    paths <- opt$search.state$ensemble.paths
    successful <- vapply(paths, `[[`, logical(1), "successful")
    errors <- vapply(paths, `[[`, character(1), "error")
    list(
        estimability=opt$search.state$estimability,
        ensemble=list(
            attempted=length(paths), successful=sum(successful),
            failed=sum(!successful),
            failure.messages=sort(table(errors[nzchar(errors)]), decreasing=TRUE)
        )
    )
}
