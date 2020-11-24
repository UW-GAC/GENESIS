varCompCI <- function(null.model, prop=TRUE){
    # check for mixed model
    if(!(null.model$model$family$mixedmodel)){
        stop("null.model is not a mixed model; this calculation can not be performed")
    }

    if(prop){
        if(null.model$model$hetResid){
            stop("Estimates of proportional variance are not supported with heterogeneous group residual variances")
        }
        ci <- matrix(NA, nrow=length(null.model$varComp), ncol=2)
        est <- null.model$varComp/sum(null.model$varComp)
        varCompCov <- null.model$varCompCov
        varCompCov[is.na(varCompCov)] <- 0
        for(i in 1:length(est)){
            deltaH <- rep(-null.model$varComp[i]/(sum(null.model$varComp)^2),length(null.model$varComp))
            deltaH[i] <- deltaH[i] + sum(null.model$varComp)/(sum(null.model$varComp)^2)
            varH <- as.numeric(crossprod(deltaH, crossprod(varCompCov, deltaH)))
            ci[i,] <- est[i] + sqrt(varH)*qnorm(c(0.025,0.975))
        }
        ci[null.model$zeroFLAG,] <- NA
        res <- as.data.frame(cbind(est, ci))
        names(res) <- c("Proportion", "Lower 95", "Upper 95")

    }else{
        ci <- null.model$varComp + sqrt(diag(null.model$varCompCov)) %o% qnorm(c(0.025,0.975))
        res <- as.data.frame(cbind(null.model$varComp, ci))
        names(res) <- c("Est", "Lower 95", "Upper 95")
    }

    print(res)
}
