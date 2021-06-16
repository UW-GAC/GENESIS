## will update this file name to nullModelInvNorm.R
## takes an existing null model (family has to be "gaussian"!), rank-normalize residuals and scale, and re-fit null model. 
## If group.idx are provided, optional scale by group.idx.
## the names of items in the list group.idx have to match the names of the corresponding variance components!

nullModelInvNorm <- function(null.model, 
                             cov.mat = NULL,
                             norm.option = c("all", "by.group"),
                             rescale = c("residSD", "none", "model"),
                             AIREML.tol = 1e-4,
                             max.iter = 100, 
                             EM.iter = 0,
                             drop.zeros = TRUE,
                             return.small = FALSE,
                             verbose = TRUE) {

    norm.option <- match.arg(norm.option)
    rescale <- match.arg(rescale)

    # checks
    if(null.model$model$family$family != "gaussian") stop("family must be gaussian")
    
    if(rescale == "model" & norm.option == "all"){
        stop("Rescaling using model variance components when rank normalizing all samples together is not implemented")
    }

    # Update null model format.
    null.model <- .updateNullModelFormat(null.model)

    resid <- null.model$fit$resid.marginal
    group.idx <- null.model$group.idx
    if(!is.null(cov.mat)){
      cov.mat <- .setCovMatNames(cov.mat)
    }

    if (norm.option == "by.group"){
        if(is.null(group.idx)) stop("Cannot rank normalize by group, missing group indices.")
        g <- length(group.idx)

        # compute variance by group
        if (rescale == "none"){
            group.vars <- rep(1, g)
        } else{
            group.vars <- rep(NA, g)
            for (i in 1:g){
                if (rescale == "model") group.vars[i] <- .averageGroupVar(null.model$varComp, cov.mat, group.idx[i])
                if (rescale == "residSD") group.vars[i] <- var(resid[group.idx[[i]]])
            }
        }

        # rank normalize and rescale
        for (i in 1:g){
            group.resids <- .rankNorm(resid[group.idx[[i]]])
            resid[group.idx[[i]]] <- group.resids*sqrt(group.vars[i])
        }
    }

    if (norm.option == "all"){

        # compute variance
        if(rescale == "none"){
            all.vars <- 1
        }else if(rescale == "residSD"){
            all.vars <- var(resid)
        }

        # rank normalize and rescale
        resid <- .rankNorm(resid)
        resid <- resid*sqrt(all.vars)
    }

    # re-fit the null model
    new.null.model <- .fitNullModel(y = resid, X = null.model$model.matrix, covMatList = cov.mat,
                                 group.idx = group.idx, family = null.model$model$family, start = null.model$varComp,
                                 AIREML.tol = AIREML.tol, max.iter = max.iter, EM.iter = EM.iter,
                                 drop.zeros = drop.zeros, return.small = return.small, verbose = verbose)

    # add any extra slots
    extra <- setdiff(names(null.model), names(new.null.model))
    new.null.model <- c(new.null.model, null.model[extra])

    # Add sample id. If the fit data frame doesn't have a sample.id column, this step does nothing.
    # In that case, it's trying to set new.null.model$fit$sample.id to NULL, so it won't exist in the
    # new fit data frame.
    new.null.model$fit$sample.id <- null.model$fit$sample.id

    # Update model strings.
    previous_model <- null.model$model$formula
    new.null.model$model$outcome <- null.model$model$outcome
    new.null.model$model$covars <- null.model$model$covars
    new.null.model$model$formula <- paste(.modelOutcomeString(null.model$model$outcome, inverse_normal=TRUE), "~",
                               strsplit(previous_model, " ~ ")[[1]][2])

    new.null.model

}


## here group.idx is a list with just one entry - a vector of the group indices.
## the name of this list entry is the name of the group, same as it appears in the varComp vector.
.averageGroupVar <- function(varComp, cov.mat = NULL, group.idx = NULL){

    if (is.null(group.idx)){
        stop("group indices are not provided, cannot calculate average variance in the group")
    }

    if (!is.list(group.idx)) stop("group.idx should be a list with one entry")
    if (length(group.idx) > 1) stop("group.idx should be a list with one entry")

    if (is.null(cov.mat)) {
        m <- 0
    } else{
        m <- length(cov.mat)
    }

    ## initialize sum of variance components
    sum.var <- 0
    if (m > 0){
        for (i in 1:m){
            sum.var <- sum.var + varComp[i]*mean(diag(cov.mat[[i]])[group.idx[[1]]])
        }
    }

    ## now add residual variance:
    sum.var <- sum.var + varComp[paste0("V_", names(group.idx)[1])]

    return(sum.var)
}



## Inverse normal transform
.rankNorm <- function(x) qnorm((rank(x) - 0.5)/length(x))
