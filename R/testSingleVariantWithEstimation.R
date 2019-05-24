# .libPaths(c("/home/students/zty/GENESISpcg", "/projects/resources/gactools/R_packages/library-3.5.1"))
# library(SeqArray)
# library(SeqVarTools)
# library(Matrix)
# library(GENESIS)
# library(devtools)
# load_all("/home/students/zty/GENESISpcg")

# set.seed(1)

calcGtildeWithW <- function(nullmod, G, r=1){
  ##Calculate Gtilde whose square is GWG
  ##Input:
  # X: design matrix
  # G: collection of single variants
  # W: diagonal matrix computed from variance components
  ##Returns:
  # Gtilde: W^(1/2)(G - X(XWX)^(-1)XWG)
  
  X <- nullmod$model.matrix
  W <- nullmod$W
  
  XWX.inv <- solve(crossprod(X,crossprod(W,X)))
  Gtilde <- G - X %*% (XWX.inv %*% crossprod(X,crossprod(W,G)))
  Gtilde <- r^(1/2)*crossprod(W^(1/2),Gtilde)
  return(Gtilde)
  }

calW <- function(nullmod){
  ###Calculate W matrix
  ##Input:
  # nullmod: null model
  ##Output:
  # W
  
  Y <- nullmod$outcome
  varComp <- nullmod$varComp
  group.idx <- nullmod$group.idx
  vmu <- nullmod$vmu
  
  m <- 1 #number of covariance matrix
  n <- length(Y)
  if (is.null(vmu)){ ## this means the family is "gaussian"
    if (is.null(group.idx)){
      diagV <- rep(varComp[m+1],n)
    } else{
      
      g <- length(group.idx)
      mylevels <- rep(NA, n)
      
      for(i in 1:g){
        mylevels[group.idx[[i]]] <- i # setting up vector of indicators; who is in which group
      }
      diagV <- (varComp[m+1:g])[mylevels]
    }
    
    W <- diag(1/diagV)
  }
  return(W)
}

estr <- function(nullmod, G,
                  rEstSize = 30){
##Estimating the variance ratio r = GPG/GWG
##Input:
  # nullmod: null model
  # G: collection of single variants
  # rEstSize: number of variants used to estimate ratio
##Returns:
  # r: Estimated ratio
  # GPG: true variance estimator
  # GWG: estimated variance estimator
  
  Y <- nullmod$outcome
  X <- nullmod$model.matrix
  W <- nullmod$W
  
  ncolG <- dim(G)[2]
  
  index <- sample(x = 1:ncolG,size = rEstSize)
  G <- G[,index] ## We choose 30 genes.
  
  Gtilde1 <- calcGtilde(nullmod, G)
  GPG <- colSums(Gtilde1^2)
  
  Gtilde2 <- calcGtildeWithW(nullmod, G)
  GWG <- colSums(Gtilde2^2)
  
  r <- mean(GPG/GWG)
  return(list(GPG,GWG,r))
}

# #####homo SKM
# file_nullmod <- "/projects/topmed/research/sparse_matrices/analysts/mconomos/HCT/pipeline_jobs/assoc_single_noInvNorm/pcrelate_all_samples_sparse5_subset_pruned_homVar/data/pcrelate_all_samples_sparse5_pruned_homVar_null_model.RData"
# nullmod <- get(load(file_nullmod))
# file_assoc <- "/projects/topmed/research/sparse_matrices/analysts/mconomos/HCT/pipeline_jobs/assoc_single_noInvNorm/pcrelate_all_samples_sparse5_subset_pruned_homVar/data/pcrelate_all_samples_sparse5_pruned_homVar_assoc_single_chr22.RData"
# assochr22 <- get(load(file_assoc))
# # file_est_r <- "/home/students/zty/src/AssociateTest/UseGWG/estimater_heter_KM.rda"
# # FileTenBlock <- "/home/students/zty/src/AssociateTest/UseGWG/GWG_homo_KM_10blocks.rda"
# FilePlotData <- "/home/students/zty/src/AssociateTest/UseGWG/GPG-GWG-homo-SKM.rda"
# FilePlot <- "/home/students/zty/src/AssociateTest/UseGWG/GPG-GWG-homo-SKM.png"
# 
# ###Prepare G matrix
# gdsfile <- "/projects/topmed/downloaded_data/IRC_freezes/freeze.5b/gds/minDP0/freeze5b.chr22.subset_fr5.pass_and_fail.gtonly.minDP0.gds"
# gds <- seqOpen(gdsfile)
# index <- get(load("/home/students/zty/src/AssociateTest/UseGWG/index.rda"))
# 
# ##build a small G matrix with more than 30 columns
# variant <- index[[2]]
# random100 <- sample(1:length(variant),size = 100)
# seqSetFilter(gds,sample.id = index[[1]],variant.id = variant[random100])
# G <- expandedAltDosage(gds, use.names=T, sparse =TRUE)
# 
# ##Choose mac > 20
# G <- G[,pmin(colSums(G), 2*nrow(G) - colSums(G))>20]
# 
# #estimate r
# W <- calW(nullmod)
# b <- estr(nullmod, G)
# 
# #estimate variance for block*N variants
# block <- 1024
# N <- 10
# k <- 0
# allGWG <- rep(0,block*N)
# 
# while(k < N){
#   seqSetFilter(gds,sample.id = index[[1]],variant.id = variant[(k*block+1):((k+1)*block)])
#   G <- expandedAltDosage(gds, use.names=T, sparse =TRUE)
#   Gtilde <- calcGtildeWithW(X, G, W)
#   GWG <- colSums(Gtilde^2)
#   allGWG[(k*block+1):((k+1)*block)] <- GWG
#   names(allGWG)[(k*block+1):((k+1)*block)] <- variant[(k*block+1):((k+1)*block)]
#   k <- k+1
# }
# 
# 
# # #############
# # ###Get p-value
# ratio <- b[[3]]
# GWG <- allGWG
# 
# Score.SE.W <- sqrt(GWG*ratio)
# Score.Stat.W <- dataGENESIS$Score[1:length(Score.SE.W)]/Score.SE.W
# Score.pval.W <- 1- pchisq((Score.Stat.W)^2, 1)
# 
# # library(ggplot2)
# # library(hexbin)
# # df <- data.frame(GPG.P = dataGENESIS$Score.pval[1:length(Score.SE.W)], GWG.P = Score.pval.W,
# #                  ratio = rs, largeratio = rsb)
# # p1 <- ggplot()+
# #   geom_hex(aes(x=-log10(GPG.P),y=-log10(GWG.P),fill = log(..count..)),data = df, bins =100)
# # ggsave(FilePlot)
# # data <- list(df,ratio)
# # save(data, file = FilePlotData)
# # # summary((-log10(df$GPG.P))/(-log10(df$GWG.P)))
# # p2 <- ggplot()+geom_histogram(aes(x=ratio),data= df)
# # ggsave("ratio-grm-small-chr22.png")
# #   geom_hex(aes(x=-log10(GPG.P),y=-log10(GWG.P),fill = log(..count..)),data = df, bins =100)
# # subindex <- sample(x = 1:dim(df)[1], size = 100)
# # subdf <- df[subindex,]
# # p3 <- ggplot()+geom_point(aes(x=-log10(GPG.P),y=-log10(GWG.P),col = largeratio),data= df)
# # ggsave("scattter-grm-small-chr22.png")
# # 
# # 
# # # p3 <- ggplot()+geom_point(aes(x=GWG,y=GPG),data= df)
# # # ggsave("GWG-GPG-grm-small-chr22.png")
# # # p2 <- ggplot()+geom_histogram(aes(x=GWG),data= df)
# # # ggsave("GWG-grm-small-chr22.png")
# # # p3 <- ggplot()+ geom_hex(aes(x=rs,y=Gs,fill = log(..count..)),data = df, bins =100)
# # # ggsave("rs-GcolSums-grm-small-chr22.png")
# # 
# # ####investigate why there are two groups
# # load('/projects/topmed/research/sparse_matrices/analysts/mconomos/HCT/data/pheno.RData')
# # GWG <- get(load("/home/students/zty/src/AssociateTest/UseGWG/GWGsmall.rda"))
# # dataGENESIS <- get(load("/home/students/zty/src/AssociateTest/20kGENESIS/data/20kGENESIS_assoc_single_chr22.RData"))
# # Score.SE.W <- sqrt(GWG*0.6588)
# # rs <- (dataGENESIS$Score.SE[1:length(Score.SE.W)])^2/GWG
# # 
# # 
# # datafull <- tabfull@data[complete.cases(tabfull@data[,5]),]
# # rownames(datafull) <- datafull$sample.id
# # table(datafull$ancestry)
# # euroid <- (datafull$sample.id)[datafull$ancestry == "European American"]
# # afrid <- datafull$sample.id[datafull$ancestry == "African American"]
# # library(SeqArray)
# # library(SeqVarTools)
# # gdsfile <- "/projects/topmed/downloaded_data/IRC_freezes/freeze.5b/gds/minDP0/freeze5b.chr22.subset_fr5.pass_and_fail.gtonly.minDP0.gds"
# # gds <- seqOpen(gdsfile)
# # index <- get(load("/home/students/zty/src/AssociateTest/UseGWG/index.rda"))
# # variant <- index[[2]]
# # block <- 1024
# # k <- 0
# # seqSetFilter(gds,sample.id = c(euroid,afrid),variant.id = variant[k*block+1:((k+1)*block)])
# # G <- expandedAltDosage(gds, use.names=T, sparse =TRUE)
# # euroG <- G[euroid,]
# # afrG <- G[afrid,]
# # euromean <- colSums(euroG)/length(euroid)
# # afrmean <- colSums(afrG)/length(afrid)
# # eurosum <- colSums(euroG)
# # afrmean <- colSums(afrG)/length(afrid)
# # dmean <- euromean - afrmean
# # 
# # library(ggplot2)
# # library(hexbin)
# # df <- data.frame(dmean = dmean, rs = rs[1:block])
# # p1 <- ggplot()+
# #   geom_hex(aes(x=dmean,y=rs,fill = log(..count..)),data = df, bins =100)
# # ggsave("meandif-rs.png")
# # 
# # meanG <- rep(0,dim(G)[2])
# # for(i in 1:dim(G)[2]){
# #   meanG[i] <- mean(G[,i])
# #   }
# # pdf("meanG-rs.pdf") 
# # # plot(dmean,rs[1:block])
# # plot(meanG,rs[1:dim(G)[2]])
# # dev.off()
# # 
# # 
# # 
# # saiges <- read.table("/home/students/zty/src/SAIGEresult_500markers.SAIGE.results.txt", header = T)
# # saige500 <- saiges$var1/saiges$var2
# # pdf("saige500-hist.pdf") 
# # hist(saige500,breaks = 30)
# # dev.off()