### new SPA function

SPA_pval <- function(score.result, nullmod, G, pval.thresh = 1){
	if (!requireNamespace("SPAtest")) stop("package 'SPAtest' must be installed to calculate SPA p-values")
	expit <- function(x){exp(x)/(1+exp(x))}

	# index for variants with Score.pval < pval.thresh
	# only bother running SPA on these variants
	idx <- which(score.result$Score.pval <= pval.thresh)

	# update columns in score.result
	setnames(score.result, "Score.pval", "SPA.pval")
        if (nrow(score.result) > 0) {
            score.result$SPA.converged <- NA
        } else {
            return(cbind(score.result, data.frame(SPA.converged=logical())))
        }

	if(length(idx) > 0 ){

		# calculate the fitted values  from the null model;
                # expit(eta) = expit(X\beta + Zb)
                if (nullmod$model$family$mixedmodel) {
                    mu <- as.vector(expit(nullmod$fit$linear.predictor))
                } else {
                    # fitted.values from glm models are already the expit.
                    mu <- as.vector(nullmod$fit$fitted.values)
                }

		# W is the diagonal of a matrix
		W <- mu*(1-mu)
		WX <- W*nullmod$model.matrix
		XWX.inv <- solve(crossprod(nullmod$model.matrix,WX))

		# original code filtered variants with MAC <= 3 here

		# loop through variants
		for(i in idx){
			# extract the genotypes
			g <- G[,i]
			# get the score
			s <- score.result$Score[i]

			# "flip" g and s if the minor allele is the reference allele
			if(mean(g) > 1){
				g <- 2-g
				s <- -s
			}
			# identify which elements in g are hom ref
			homRefSet <- which(g == 0)

			# compute adjusted genotype values
			# G - X(X'WX)^{-1}(X'WG)
			g <- as.vector(g - tcrossprod(nullmod$model.matrix, crossprod(crossprod(WX, g), XWX.inv)))

			# compute variance ratio (similar to SAIGE estimator)
			GPG <- score.result$Score.SE[i]^2
			GWG <- sum(W*g^2)
			r <- GPG/GWG

			# get inputs to SPAtest function
			g <- sqrt(r)*g
			qtilde <- as.numeric(s + crossprod(g, mu))

			# compute SPA p-value
			if(length(homRefSet)/length(g) < 0.5){
	        	tmp <- SPAtest::Saddle_Prob(q = qtilde, mu = mu, g = g, alpha=5e-8, output = "P")
	        }else{
	        	tmp <- SPAtest::Saddle_Prob_fast(	q = qtilde, mu = mu, g = g, alpha = 5e-8, output = "P",
	        										gNA = g[homRefSet], gNB = g[-homRefSet],
	        										muNA = mu[homRefSet], muNB = mu[-homRefSet])
	        }

	        # add in results
	        score.result$SPA.converged[i] <- tmp$Is.converge
	        if(tmp$Is.converge){
	        	score.result$SPA.pval[i] <- tmp$p.value
	        }
		}
	}

	return(score.result)
}

### some code for vectorized version of for loop above
# Gadj <- G[,idx] - tcrossprod(nullmod$model.matrix, crossprod(crossprod(WX, G[,idx]), XWX.inv))

# # compute variance ratio (similar to SAIGE estimator)
# GPG <- score.result$Score.SE[idx]^2
# GWG <- colSums(W*Gadj^2)
# r <- GPG/GWG

# # compute inputs to SPAtest function
# Gadj <- t(t(Gadj)*sqrt(r))
# qtilde <- score.result$Score[idx] + crossprod(Gadj, mu)
