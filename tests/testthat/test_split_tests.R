context("split assoc tests")
library(SeqVarTools)

test_that("genIDList - check errors", {
  dat <- pData(sampleData(.testData()))
  id.var="sample.id"
  group.var="sex"
  mat <- as.matrix(cbind(dat$sample.id, dat$sex))
  colnames(mat) <- c("sample.id", "sex")
  expect_error(genIDList(mat, id.var=id.var, group.var=group.var))
  expect_error(genIDList(dat, id.var="sample", group.var=group.var))
  expect_error(genIDList(dat, id.var=id.var, group.var="bmi"))
})


test_that("genIDList", {
  dat <- pData(sampleData(.testData()))
  id.var="sample.id"
  n_cat <- length(unique(dat$sex))
  test_list <- genIDList(dat, id.var=id.var, group.var="sex")
  expect_equal(length(test_list), n_cat)
  i=1
  current_union=test_list[[i]]
  current_sum=length(test_list[[i]])
  while (i < length(test_list)){
    current_union=union(current_union, test_list[[i+1]])
    current_sum = current_sum + length(test_list[[i+1]])
    i=i+1
  }
  expect_equal(length(current_union), current_sum)
})

###nullModelSplit currently assumes null model is a mixed model - should put logical on nullmod$family to see if resid.conditional should be subsetted
test_that("nullModelSplit", {
  test_nullmod <- .testNullmod(n=100, MM=TRUE, binary=TRUE)
  test_nullmod$sample.id <- paste0("A", as.character(seq(1,100)))
  id_list = list()
  id_list[['g1']] <- test_nullmod$sample.id[1:40]
  id_list[['g2']] <- test_nullmod$sample.id[41:100]
  split_nullmod_all = nullModelSplit(test_nullmod, id_list, keep.all=TRUE)
  split_nullmod_groups= nullModelSplit(test_nullmod, id_list, keep.all=FALSE)
  check_values <- c('outcome', 'workingY', 'resid.conditional', 'fitted.values')
  for (g in names(id_list)){
    expect_true(all(vapply(check_values, function(x, data) length(data[[x]]) == length(data[['sample.id']]), logical(1), data=split_nullmod_all[[g]])))
  }
})


test_that("assocTestSingleSplit", {
  svd <- .testData()
  iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
  nullmod <- fitNullModel(iterator, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
  sample_ids <- nullmod$sample.id
  id_list = list()
  id_list[['g1']] <- sample_ids[1:40]
  id_list[['g2']] <- sample_ids[41:length(sample_ids)]
  split_nullmod = nullModelSplit(nullmod, id_list, keep.all=FALSE)
  assoc = assocTestSingle(iterator, split_nullmod[['g2']], test="BinomiRare", verbose=FALSE)
  resetIterator(iterator, verbose=FALSE)
  split_assoc = assocTestSingleSplit(iterator, nullmod, id.list=id_list, test="BinomiRare", verbose=FALSE)
  expect_equal(assoc, split_assoc[['g2']])
  seqClose(svd)
  }
)


test_that("assocTestAggregateSplit", {
  svd <- .testData()
  gr <- granges(svd)
  vi <- variantInfo(svd, alleles=TRUE)
  mcols(gr) <- vi[,c("ref", "alt")]
  grl <- GRangesList(gr[c(1,1,2)])
  iterator <- SeqVarListIterator(svd, variantRanges=grl, verbose=FALSE)
  
  nullmod <- fitNullModel(iterator, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
  sample_ids <- nullmod$sample.id
  id_list = list()
  id_list[['g1']] <- sample_ids[1:40]
  id_list[['g2']] <- sample_ids[41:length(sample_ids)]
  split_nullmod = nullModelSplit(nullmod, id_list, keep.all=FALSE)
  
  assoc <- assocTestAggregate(iterator,split_nullmod[['g2']], test="BinomiRare", verbose=FALSE)
  resetIterator(iterator, verbose=FALSE)
  split_assoc <- assocTestAggregateSplit(iterator, nullmod, id.list=id_list, test="BinomiRare", keep.all=FALSE, verbose=FALSE)
  expect_equal(assoc, split_assoc[['g2']])
  seqClose(svd)
}
)


test_that("assocTestAggregateSplit - CMP", {
  svd <- .testData()
  gr <- granges(svd)
  vi <- variantInfo(svd, alleles=TRUE)
  mcols(gr) <- vi[,c("ref", "alt")]
  grl <- GRangesList(gr[c(1,1,2)])
  iterator <- SeqVarListIterator(svd, variantRanges=grl, verbose=FALSE)
  
  nullmod <- fitNullModel(iterator, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
  sample_ids <- nullmod$sample.id
  id_list = list()
  id_list[['g1']] <- sample_ids[1:40]
  id_list[['g2']] <- sample_ids[41:length(sample_ids)]
  split_nullmod = nullModelSplit(nullmod, id_list, keep.all=FALSE)
  
  assoc <- assocTestAggregate(iterator,split_nullmod[['g2']], test="CMP", verbose=FALSE)
  resetIterator(iterator, verbose=FALSE)
  split_assoc <- assocTestAggregateSplit(iterator, nullmod, id.list=id_list, test="CMP", keep.all=FALSE, verbose=FALSE)
  expect_equal(assoc, split_assoc[['g2']])
  seqClose(svd)
  
}
)

test_that("assocTestAggregateSplit - window", {
  svd <- .testData()
  seqSetFilterChrom(svd, include=1, verbose=FALSE)
  iterator <- SeqVarWindowIterator(svd, windowSize=5e5, windowShift=2.5e5, verbose=FALSE)
  nullmod <- fitNullModel(iterator, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
  
  nullmod <- fitNullModel(iterator, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
  sample_ids <- nullmod$sample.id
  id_list = list()
  id_list[['g1']] <- sample_ids[1:40]
  id_list[['g2']] <- sample_ids[41:length(sample_ids)]
  split_nullmod = nullModelSplit(nullmod, id_list, keep.all=FALSE)
  
  assoc <- assocTestAggregate(iterator,split_nullmod[['g2']], test="CMP", verbose=FALSE)
  resetIterator(iterator, verbose=FALSE)
  split_assoc <- assocTestAggregateSplit(iterator, nullmod, id.list=id_list, test="CMP", keep.all=FALSE, verbose=FALSE)
  expect_equal(assoc, split_assoc[['g2']])
  seqClose(svd)
})

