context("test split utils")
library(SeqVarTools)


test_that('detect variant ID variable', {
  svd <- .testData()
  iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
  nullmod <- fitNullModel(iterator, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
  assoc <- assocTestSingle(iterator,nullmod, test="BinomiRare", verbose=FALSE)
  expect_false(.detectVarNames(assoc, "variantID"))
  expect_true(.detectVarNames(assoc, "variant.id"))
  expect_true(.detectVarNames(assoc, c("variantID", 'variant.id', "abc")))
  expect_false(.detectVarNames(assoc, "abc"))
  seqClose(svd)
}
)



test_that('mergeNullModelBR', {
  svd <- .testData()
  grm <- .testGRM(svd)
  sampleData <- sampleData(svd)
  nullmod.list <- vector(mode="list", length=2)#fitNullModel(svd, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
  sample_ids <- sampleData$sample.id
  id_list = list()
  id_list[['g1']] <- sample_ids[1:40]
  id_list[['g2']] <- sample_ids[41:length(sample_ids)]
  nullmod.list[[1]] <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), cov.mat=grm[id_list[['g1']], id_list[['g1']]], family="binomial", verbose=FALSE)
  nullmod.list[[2]] <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), cov.mat=grm[id_list[['g2']], id_list[['g2']]], family="binomial", verbose=FALSE)
  merged_null <- mergeNullModelBR(nullmod.list, svd)
  
  expect_equivalent(sampleData$status, merged_null$outcome)
  expect_equal(sampleData$sample.id, merged_null$sample.id)
})


test_that('check matchSignifHits rows, length of findIncompleteHits, and runSplitSubset rows', {
  svd <- .testData()
  iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
  nullmod <- fitNullModel(iterator, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
  sample_ids <- nullmod$sample.id
  id_list = list()
  id_list[['g1']] <- sample_ids[1:40]
  id_list[['g2']] <- sample_ids[41:length(sample_ids)]
  split_assoc = assocTestSingleSplit(iterator, nullmod, id.list=id_list, test="BinomiRare", verbose=FALSE)
  matched_assoc_list <- matchSignifHits(res.list=split_assoc, threshold=0.1)
  matched_assoc_df <- matchSignifHits(split_assoc, threshold=0.1, return.df=TRUE)
  ##expect same number of rows for results when returning list or dataframe
  expect_equal(nrow(dplyr::bind_rows(matched_assoc_list)), nrow(matched_assoc_df))
  expect_equal(length(unique(matched_assoc_list[['g1']]$variant.id)), sum(matched_assoc_list[['g1']]$ref.group=='all')) ##variants in each group result should have a corresponding result in 'all'
  incomplete_hits <- findIncompleteHits(matched_assoc_list)
  expect_equal(length(incomplete_hits), length(findIncompleteHits(matched_assoc_df)))
  resetIterator(iterator, verbose=FALSE)
  recreateIterator(iterator, annot=sampleData(svd), incomplete.variants=incomplete_hits, verbose=FALSE)
  expect_equal(length(seqGetData(iterator, "variant.id")), length(incomplete_hits))
  seqClose(svd)
})


test_that('runSplitSubset',{
  svd <- .testData()
  variants <- seqGetData(svd, "variant.id")[1:30]
  seqSetFilter(svd, variant.id=variants, verbose=FALSE)
  nullmod <- fitNullModel(sampleData(svd), outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
  sample_ids <- nullmod$sample.id
  id_list = list()
  id_list[['g1']] <- sample_ids[1:40]
  id_list[['g2']] <- sample_ids[41:length(sample_ids)]
  iterator <- SeqVarBlockIterator(svd, verbose=FALSE)
  assoc_subset <- runSplitSubset(iterator, nullmod, id.list=id_list, test="BinomiRare", verbose=FALSE)
  expect_equal(length(unique(unlist(lapply(assoc_subset, nrow)))), 1)
  seqClose(svd)
}
)




