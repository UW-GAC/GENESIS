

# nullmod: null model results 
# obj.glm.null: glm object, for example obj.glm.null = glm( y~x1+x2 , data=data, family=binomial)
# score.result is the result from score test
# SNP is the name of the markers
# N is the sample size
SAIGE_Pvalue <- function(nullmod, score.result, geno, SNP=NULL, N=NULL){
#    if (!requireNamespace("SPAtest")) stop("package 'SPAtest' must be installed to calculate SAIGE p-values") ##already in pkg NAMESPACE
    if(is.null(N[1])){N <- rep(nrow(geno),ncol(geno))}
    # calculate allele frequency
    AF <- colMeans(geno)/2
    obj.noK <- SPAtest::ScoreTest_wSaddleApprox_NULL_Model( nullmod$outcome~nullmod$model.matrix)
    y <- obj.noK$y
    expit <- function(x){exp(x)/(1+exp(x))}
    mu <- expit(nullmod$fitted.values)
    
    # calculate SAIGE pvalues by for loop
    out <- NULL
    for (i in 1:ncol(geno)){
        #obtain a vector with genotypes and a list with association results for the ith marker
        G <- geno[,i]
        score.resultList <- list(SNP=SNP[i], N=N[i], AF=AF[i],
                                 SCORE=score.result$Score[i], VAR=score.result$Score.SE[i]^2, PVAL=score.result$Score.pval[i])
        
        if(score.resultList$AF == 0 | score.resultList$AF == 1 | (2*score.resultList$AF*score.resultList$N) <= 3 ){
            result <- c(score.resultList, list(PVAL.saige=NA, PVAL.NA=NA, SPA.Is.converge=NA))
        }else{
            #impute the missing genotype
            G[which(is.na(G))] <- 2*(score.resultList$AF[i])
            G <- as.numeric(G)
            result <- getSAIGEPvalfromGMMATOutput(y,mu,score.resultList,G,obj.noK)
        }
        out <- rbind(out, unlist(result))
    }
    return(data.frame(out)) 
}



#This function applies saddlepoint approximation given the GMMAT results. Tstat is SCORE in GMMAT output. var1 is VAR in GMMAT output
scoreTest_SPAGMMAT_binaryTrait=function(g, AC, NAset, y, mu, Tstat, var1, Cutoff = 2){
    q <- crossprod(g, y)
    m1 <- crossprod(g, mu)
    var2 <- crossprod(mu*(1-mu), g*g)
    qtilde <- (Tstat/sqrt(AC))/sqrt(var1/AC) * sqrt(var2) + m1
    if(length(NAset)/length(g) < 0.5){
        out1 <- SPAtest:::Saddle_Prob(q=qtilde, mu = mu, g = g, Cutoff = Cutoff, alpha=5*10^-8, nodes.fixed=NULL, nodes.init=NULL)
    }else{
        out1 <- SPAtest:::Saddle_Prob_fast(q=qtilde,g = g, mu = mu, gNA = g[NAset], gNB = g[-NAset], muNA = mu[NAset], muNB = mu[-NAset], Cutoff = Cutoff, alpha = 5*10^-8, output="P", nodes.fixed=NULL, nodes.init=NULL)
    }
    
    out <- list(PVAL.saige=out1$p.value, PVAL.NA=out1$p.value.NA, SPA.Is.converge=out1$Is.converge)
    return(out) 
}



##SNP     N       AF      SCORE   VAR     PVAL

#This function obtain the SAIGE/SPA p value for each genetic marker
#G: a vector containing raw gentoypes for the marker
#GMMATAssocResultList: a list containing association results for the marker 
getSAIGEPvalfromGMMATOutput <- function(y, mu, GMMATAssocResultList, G, obj.noK){
    
    mu.a<-as.vector(mu)
    
    if(GMMATAssocResultList$PVAL <= 0.05 & GMMATAssocResultList$AF*2*GMMATAssocResultList$N > 2){ 
        G0 <- matrix(G, ncol = 1)
        N <- GMMATAssocResultList$N
        if(GMMATAssocResultList$AF > 0.5){
            G0 <- 2-G0
            AC2 <- 2*N - (GMMATAssocResultList$AF*2*N)
        }else{
            AC2 <- GMMATAssocResultList$AF*2*N
        }
        XVG0 <- as.matrix(obj.noK$XV)%*%matrix(G0,ncol=1)
        G1 <- G0  -  (obj.noK$XXVX_inv)%*%XVG0 # G is X adjusted
        g <- G1/sqrt(AC2)
        NAset <- which(G0==0)
        pvalue.SPA <- scoreTest_SPAGMMAT_binaryTrait(g, AC2, NAset, y, mu.a, GMMATAssocResultList$SCORE,GMMATAssocResultList$VAR)
    }else{
        pvalue.SPA <- list(PVAL.saige=NA, PVAL.NA=NA, SPA.Is.converge=NA)	
    }
    return(c(GMMATAssocResultList,pvalue.SPA))
}
