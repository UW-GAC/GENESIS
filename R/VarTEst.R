#### This code is for estimating the vairance of score test statistic 
#### Based on SAGIE's method, created by Tinayu.

var_ratio_est <- function(nullmod,covMatList){
  mu <- nullmod$fitted.values
  hattau <- nullmod$varComp[1]
  psi <- covMatList[[1]]
  barw <- mean((mu*(1-mu))^(-1))
  barw <- nullmod$varComp[2]
  lambda <- eigen(mypcrel$kinship)$value
  r1 <- sum(lambda/(barw+hattau*lambda))
  r2 <- sum(lambda)/barw
  hatr <- r1/r2
  return(hatr)
}

#testGeno.R
#start with the first one
#.testGenoSingleVarScore #gives test statistic
#res is the final result .
##speed and result comparison
#replace calcGtilde with SAIGE Gtilde
## fstSKAT for fat eigenvale
##finding structure with randomness: probebilistic algorithm for...
#by Halko