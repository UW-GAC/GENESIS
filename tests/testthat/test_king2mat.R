context("king2mat tests")

test_that("king2mat", {
	kinfile <- system.file("extdata", "MXL_ASW.kin", package="GENESIS")
	kin0file <- system.file("extdata", "MXL_ASW.kin0", package="GENESIS")

	expect_error(king2mat(kin0file, kinfile, iids = 1:100), "Some of the provided iids are not in the KING output")
})
