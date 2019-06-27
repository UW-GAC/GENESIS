calcGtildeWithW <- function(nullmod, G, r=1){
  ## Calculate Gtilde whose square is GWG
  ## replace C = sqrt(Sigma^{-1}) with W^{1/2}

  X <- nullmod$model.matrix
  W <- nullmod$W
  
  # XWX.inv <- solve(crossprod(X,crossprod(W,X)))
  # Gtilde <- G - X %*% (XWX.inv %*% crossprod(X,crossprod(W,G)))
  # Gtilde <- r^(1/2)*crossprod(W^(1/2),Gtilde)

  # W is the diagonal of a matrix
  WX <- W*X
  XWX.inv <- solve(crossprod(X,WX))
  # G - X(X'WX)^{-1}(X'WG)
  Gtilde <- G - tcrossprod(X, crossprod(crossprod(WX, G), XWX.inv))
  Gtilde <- (sqrt(r)*sqrt(W))*Gtilde
  return(Gtilde)
}


# ## Matt - don't need this function because we save W in nullmod from .computeSigmaQuantities 6/27/2019
# calW <- function(nullmod){
#   ###Calculate W matrix

#   Y <- nullmod$outcome
#   varComp <- nullmod$varComp
#   group.idx <- nullmod$group.idx
#   vmu <- nullmod$vmu
  
#   m <- 1 #number of covariance matrix
#   n <- length(Y)
#   if (is.null(vmu)){ ## this means the family is "gaussian"
#     if (is.null(group.idx)){
#       diagV <- rep(varComp[m+1],n)
#     } else{
      
#       g <- length(group.idx)
#       mylevels <- rep(NA, n)
      
#       for(i in 1:g){
#         mylevels[group.idx[[i]]] <- i # setting up vector of indicators; who is in which group
#       }
#       diagV <- (varComp[m+1:g])[mylevels]
#     }
    
#     W <- 1/diagV
#   }
#   return(W)
# }

estVarRatio <- function(gdsobj, null.model, nvar = 100, MAC.thresh = 20,
                        sparse = TRUE, imputed = FALSE, 
                        verbose = TRUE){

  # don't use sparse matrices for imputed dosages
  if (imputed) sparse <- FALSE

  # coerce null.model if necessary
  if (sparse) null.model <- .nullModelAsMatrix(null.model)

  # filter samples to match null model
  sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose)

  # variant IDs
  varid <- seqGetData(gdsobj, 'variant.id')

  # sample variants iteratively until all meet MAC.thresh
  var.keep <- NULL
  var.drop <- NULL
  geno <- NULL
  while(length(var.keep) != nvar){
    var.try <- sample(varid[!(varid %in% c(var.keep, var.drop))], size = nvar - length(var.keep))
    seqSetFilter(gdsobj, variant.id = var.try)

    # read in genotypes
    if (!imputed) {
        geno.try <- expandedAltDosage(gdsobj, use.names=TRUE, sparse=sparse)[sample.index,,drop=FALSE]
    } else {
        geno.try <- imputedDosage(gdsobj, use.names=TRUE)[sample.index,,drop=FALSE]
    }
    # compute MAC
    alt.cnt <- colSums(geno.try)
    n.obs <- colSums(!is.na(geno.try))
    mac <- round(pmin(alt.cnt, 2*n.obs - alt.cnt))

    # update
    var.keep <- append(var.keep, var.try[mac >= MAC.thresh])
    var.drop <- append(var.drop, var.try[mac < MAC.thresh])
    geno <- cbind(geno, geno.try[, mac >= MAC.thresh, drop = FALSE])
  }

  # compute MAC for final set
  alt.cnt <- colSums(geno)
  n.obs <- colSums(!is.na(geno))
  mac <- round(pmin(alt.cnt, 2*n.obs - alt.cnt))

  # mean impute missing values
  if(any(n.obs < nrow(geno))){
    freq <- 0.5*colMean(geno)
    geno <- .meanImpute(geno, freq)
  }

  # estimate ratio
  est <- .estR(nullmod = null.model, G = geno)

  # table of underlying data
  tab <- data.frame(variant.id = colnames(geno),
                    n.obs = n.obs,
                    MAC = mac,
                    r.variant = est$r.variant,
                    GPG = est$GPG,
                    GWG = est$GWG)

  return(list(r = est$r, cv = est$cv, tab = tab))
}

.estR <- function(nullmod, G){
  ##Estimating the variance ratio r = GPG/GWG
  Gtilde1 <- calcGtilde(nullmod, G)
  GPG <- colSums(Gtilde1^2)
  
  Gtilde2 <- calcGtildeWithW(nullmod, G)
  GWG <- colSums(Gtilde2^2)
  
  r.variant <- GPG/GWG
  r <- mean(r.variant)
  cv <- sd(r.variant)/r
  return(list(r = r, cv = cv, r.variant = r.variant, GPG = GPG, GWG = GWG))
}



gdsobj <- seqOpen('/projects/topmed/qc/freeze.6/gds_subset/ld_pruned_maf01.gds')
null.model <- getobj('/projects/topmed/research/sparse_matrices/compare_methods/nullmod/hct/chol_dkm_het/data/chol_dkm_het_null_model.RData')

out <- estVarRatio(gdsobj = gdsobj, null.model = null.model)

