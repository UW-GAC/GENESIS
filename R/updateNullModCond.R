
### fit a conditional model based on a previously fit model. Default does not re-estimate variance components, but it may.
updateNullModCond <- function(nullmod, G, covMatList = NULL,  AIREML.tol = 1e-6, max.iter = 100, drop.zeros = TRUE, verbose = TRUE){

    ## a few checks that may be transfered to wrapper function:
    if (nullmod$hetResid & is.null(nullmod$group.idx)) stop("group indices are required for updating the null model")

    #if (updateVarComp){ ## this check may be pulled out for wrapper function.
    if (is.null(covMatList)) stop("covMatList is needed for udpating variance components")
    if (nullmod$hetResid & is.null(nullmod$group.idx)) stop("group indices are required for updating variance components")
    #}

    if (is.null(colnames(G))) colnames(G) <- paste0("var_", 1:ncol(G))
    X = cbind(nullmod$model.matrix, G)

    if (!nullmod$family$mixedmodel){ ## if it is not a mixed model, re-fit the model. (This includes heterogeneous residuals).

        return(.fitNullModel(nullmod$fit$outcome, X, covMatList = NULL,
                             group.idx = nullmod$group.idx, family = nullmod$family$family))
    }


    ### if the function reached this far, nullmod is a mixed model!
    ### Re-fit the model with the new design matrix and start point the varComp
    ## from the provided nullmod object.

    new.nullmod <- .fitNullModel(nullmod$fit$outcome, X, covMatList = covMatList,
                                 group.idx = nullmod$group.idx, family = nullmod$family$family, start = nullmod$varComp,
                                 AIREML.tol = AIREML.tol, max.iter= max.iter, drop.zeros = drop.zeros, verbose = verbose)

    ## add any extra slots
    extra <- setdiff(names(nullmod), names(new.nullmod))
    new.nullmod <- c(new.nullmod, nullmod[extra])
    # Update the model string.
    # This should eventually be done in the wrapper for this function, once it exists.
    previous_covar_string <- .modelCovarString(new.nullmod$covars)
    new.nullmod$covars <- c(new.nullmod$covars, colnames(G))
    new_covar_string <- .modelCovarString(new.nullmod$covars)
    new.nullmod$model <- sub(previous_covar_string, new_covar_string, new.nullmod$model)
    return(new.nullmod)
}
