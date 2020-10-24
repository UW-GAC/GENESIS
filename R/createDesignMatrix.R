
createDesignMatrix <- function(x, outcome, covars=NULL, group.var=NULL) {

    # remake factors to remove any missing levels
    for (f in union(covars, group.var)) {
        if(is.factor(x[[f]]) & !is.ordered(x[[f]])) {
            x[[f]] <- factor(x[[f]])
        }
    }
    
    if (!is.null(covars)) {
        model.formula <- as.formula(paste(outcome, "~", paste(covars, collapse="+")))
        # allow interactions
        covars <- unique(unlist(strsplit(covars,"[*:]")))
    } else {
        model.formula <- as.formula(paste(outcome, "~", 1))
    }
    x <- x[, unique(c(outcome, covars, group.var)), drop=FALSE]
    ##    x <- x[complete.cases(x),,drop=FALSE]
    ## the above command fails if there is more than one trait
    ## copied the command above to be executed below in the case of only one trait
    
    # group index
    if (!is.null(group.var)) {
        group.idx <- .indexList(x[[group.var]])
    } else {
        group.idx <- NULL
    }
    
    # outcome vector - preserve column name
    #y <- x[[outcome]]
    ## don't want this command if there is more than one trait
    ## i think the below if/else is enough and we don't need this command, even if only one trait
    ## need to double check though
    
    # make y one long vector if there are multiple outcomes
    # could optimize this for large sample sizes :)
    if(length(outcome)>1){
        y <- x[,outcome]
        outc <- length(outcome)
        newY <- y[,1]
        for(i in 2:length(outcome)){
            newY <- c(newY,y[,i])
        }
        y <- newY
    }else{
        x <- x[complete.cases(x),,drop=FALSE]
        y <- as.matrix(x[,outcome,drop=FALSE])
    }
    
    # create design matrix    
    X <- model.matrix(model.formula, data=x)
    # check for columns of all the same value (except the intercept)
    dropcol <- append(FALSE, apply(X[,-1,drop=FALSE], 2, var) == 0)
    if (sum(dropcol) > 0) {
        message("Covariates ",paste(colnames(X)[dropcol], collapse = ", "), " have only 1 value: they have been removed from the model")
        X <- X[,!dropcol,drop=FALSE]
    }

    list(y=y, X=X, group.idx=group.idx)
}

.indexList <- function(x) {
    groups <- unique(x)
    lapply(setNames(groups, groups), function(g) which(x == g))
}
