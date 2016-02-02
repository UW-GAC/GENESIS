varCompCI <- function(nullMMobj, prop=TRUE){
    if(prop){
        if(nullMMobj$hetResid){ 
            stop("Estimates of proportional variance are not supported with heterogeneous group residual variances")
        }
        ci <- matrix(NA, nrow=length(nullMMobj$varComp), ncol=2)
        est <- nullMMobj$varComp/sum(nullMMobj$varComp)
        varCompCov <- nullMMobj$varCompCov
        varCompCov[is.na(varCompCov)] <- 0
        for(i in 1:length(est)){
            deltaH <- rep(-nullMMobj$varComp[i]/(sum(nullMMobj$varComp)^2),length(nullMMobj$varComp))
            deltaH[i] <- deltaH[i] + sum(nullMMobj$varComp)/(sum(nullMMobj$varComp)^2)
            varH <- crossprod(deltaH, crossprod(varCompCov, deltaH))
            ci[i,] <- est[i] + sqrt(varH)*qnorm(c(0.025,0.975))
        }
        ci[nullMMobj$zeroFLAG,] <- NA
        res <- as.data.frame(cbind(est, ci))
        names(res) <- c("Proportion", "Lower 95", "Upper 95")
        
    }else{
        ci <- nullMMobj$varComp + sqrt(diag(nullMMobj$varCompCov)) %o% qnorm(c(0.025,0.975))
        res <- as.data.frame(cbind(nullMMobj$varComp, ci))
        names(res) <- c("Est", "Lower 95", "Upper 95")
    }
    
    print(res)
}
