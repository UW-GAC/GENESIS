test_king2mat <- function(){
	kinfile <- system.file("extdata", "MXL_ASW.kin", package="GENESIS")
	kin0file <- system.file("extdata", "MXL_ASW.kin0", package="GENESIS")

	obs <- tryCatch( king2mat(kin0file, kinfile, iids = 1:100), error=conditionMessage)
	checkIdentical(obs, "Some of the provided iids are not in the KING output")
}