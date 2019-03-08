
### takes an existing linear model (family has to be "gaussian"!), rank-normalize residuals and scale, and re-fit null model. If group.idx are provided, optional scale by group.idx. 

## the names of items in the list group.idx have to match the names of the corresponding variance components!

updateNullModOutcome <- function(nullmod, 
                                covMatList = NULL, 
                                rankNorm.option = c("by.group", "all"), 
                                rescale = c("none", "model", "residSD"), 
                                AIREML.tol = 1e-6, 
                                max.iter = 100, 
                                verbose = TRUE){

    rankNorm.option <- match.arg(rankNorm.option)
    rescale <- match.arg(rescale)
    
    if (nullmod$family$family != "gaussian") stop("Family must be gaussian")
    
    resid <- nullmod$resid.marginal
    
    group.idx <- nullmod$group.idx
    if (!is.null(group.idx)) {
        g <- length(group.idx)
    } else {
        g <- 1
    }

    if(!is.null(covMatList)){
        if (!is.list(covMatList)){
            covMatList <- list(A = covMatList)
        }
    }
    
    ## checks that may be put into wrapper:
    # if ((rescale != "none") & is.null(group.idx)) stop("Rescaling is only done by groups, and group indices are missing.") 
    if(rescale == "model" & rankNorm.option == "all") stop("Rescaling using model variance components when rank normalizing all samples together is not implemented")
    
    if (rankNorm.option == "by.group"){
        if(is.null(group.idx)) stop("Cannot rank normalize by group, missing group indices.")
        
        # compute variance by group
        if (rescale == "none"){
            group.vars <- rep(1, g)
        } else{
            group.vars <- rep(NA, g)
            for (i in 1:g){
                if (rescale == "model") group.vars[i] <- .averageGroupVar(nullmod$varComp, covMatList, group.idx[i])
                if (rescale == "residSD") group.vars[i] <- var(resid[group.idx[[i]]])
            }
        }

        # rank normalize and rescale
        for (i in 1:g){
            group.resids <- .rankNorm(resid[group.idx[[i]]])
            resid[group.idx[[i]]] <- group.resids*sqrt(group.vars[i])    
        }
    }
    
    if (rankNorm.option == "all"){

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

    # old code for rankNorm.option == "all"
    # for (i in 1:g){
    #     resid[group.idx[[i]]] <- resid[group.idx[[i]]]*sqrt(group.vars[i])/sd(resid[group.idx[[i]]])
    # }
    
    ### now re-fit the null model:
    
    new.nullmod <- .fitNullModel(y = resid, X = nullmod$model.matrix, covMatList = covMatList,
                                 group.idx = group.idx, family = "gaussian", start = nullmod$varComp, 
                                 AIREML.tol = AIREML.tol, max.iter = max.iter, drop.zeros = TRUE, 
                                 verbose = verbose)

    ## add any extra slots
    extra <- setdiff(names(nullmod), names(new.nullmod))
    new.nullmod <- c(new.nullmod, nullmod[extra])
    return(new.nullmod)
}


## here group.idx is a list with just one entry - a vector of the group indices.
## the name of this list entry is the name of the group, same as it appears in the varComp vector. 
.averageGroupVar <- function(varComp, covMatList = NULL, group.idx = NULL){

    if (is.null(group.idx)){
        stop("group indices are not provided, cannot calculate average variance in the group")
    } 
    
    if (!is.list(group.idx)) stop("group.idx should be a list with one entry")
    if (length(group.idx) > 1) stop("group.idx should be a list with one entry")
    
    if (is.null(covMatList)) {
        m <- 0	
    } else{
        m <- length(covMatList)
    }
    
    ## initialize sum of variance components
    sum.var <- 0
    if (m > 0){
        for (i in 1:m){
            sum.var <- sum.var + varComp[i]*mean(diag(covMatList[[i]])[group.idx[[1]]])
        }
    }
    
    ## now add residual variance:
    sum.var <- sum.var + varComp[paste0("V_", names(group.idx)[1])]
    
    return(sum.var)
}



## Inverse normal transform
.rankNorm <- function(x) qnorm((rank(x) - 0.5)/length(x))
