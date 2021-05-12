
test_that("variant specific inflation factor functions",{
  # make sure errors are checked
  freq_vec <- c(0.1, 0.2)
  ns <- c(2000, 5000, 100)
  sigma_sqs <- c(1, 1, 2)
  expect_error(computeVSIF(freq = freq_vec, ns, sigma_sqs))
  
  
  freq_vec <- c(0.1, 0.2, 0.5)
  names(freq_vec) <- names(ns) <- c("g1", "g2", "g3")  
  names(sigma_sqs) <-c("g1", "g2", "g4")  
  expect_error(computeVSIF(freq = freq_vec, ns, sigma_sqs))
  
  names(sigma_sqs) <-c("g1", "g2", "g3")  
  res <- computeVSIF(freq = freq_vec, ns, sigma_sqs)
  res[1,] <- round(res[1,],2)
  expect_equal(res, data.frame(SE_true = 1.9, SE_naive = 1.89, Inflation_factor = 1.01))
  
  freq_mat <- matrix(c(0.1, 0.2, 0.5, 0.1, 0.01, 0.5), nrow = 2, byrow = TRUE)
  colnames(freq_mat) <- names(sigma_sqs)
  res2 <- computeVSIF(freq = freq_mat, ns, sigma_sqs)
  expect_true(res2$Inflation_factor[2] > res2$Inflation_factor[1])
  
})


test_that("variant specific inflation factors using null model", {
  n <- 1000
  set.seed(22); 
  outcome <- c(rnorm(n*0.28, sd =1), rnorm(n*0.7, sd = 1), rnorm(n*0.02, sd = sqrt(2)) )
  dat <- data.frame(sample.id=paste0("ID_", 1:n),
                    outcome = outcome,
                    b=c(rep("g1", n*0.28), rep("g2", n*0.7), rep("g3", n*0.02)),
                    stringsAsFactors=FALSE)
  dat <- Biobase::AnnotatedDataFrame(dat)
  nm <- fitNullModel(dat, outcome="outcome", covars="b", verbose=FALSE)
  freq_vec <- c(0.1, 0.2, 0.5)
  names(freq_vec) <- c("g1", "g2", "g3")
  group_var_vec <- dat$b
  names(group_var_vec) <- dat$sample.id
  
  
  res <- computeVSIFNullModel(nm, freq_vec, group_var_vec)
  expect_true(res$Inflation_factor > 1)
  
  freq_mat <- matrix(c(0.1, 0.2, 0.5, 0.1, 0.01, 0.5), nrow = 2, byrow = TRUE)
  colnames(freq_mat) <- c("g1", "g2", "g3")
  res2 <- computeVSIFNullModel(nm, freq_mat, group_var_vec)
  expect_true(nrow(res2) == 2)
  
  # now check errors
  expect_error(computeVSIFNullModel(nm, freq_mat[,1:2], group_var_vec))
  
  freq_mat2 <- freq_mat
  colnames(freq_mat2)[1] <- "gf"
  expect_error(computeVSIFNullModel(nm, freq_mat2, group_var_vec))
    
  # some people not available in group_var_vec:
  expect_error(computeVSIFNullModel(nm, freq_mat2, group_var_vec[1:100]))  
  
  # additional group level in group_var_vec:
  group_var_vec2 <- group_var_vec
  group_var_vec2[100] <- "gf"
  expect_error(computeVSIFNullModel(nm, freq_mat2, group_var_vec2)) 
  
  # a person in group_var_vec that did not participate in the null model-- 
  # function should work
  group_var_vec2 <- c(group_var_vec, c(ID_another = "g3"))
  res <- computeVSIFNullModel(nm, freq_mat, group_var_vec2)
  expect_true(nrow(res) == 2)
  
})
