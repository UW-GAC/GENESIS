king2mat <-
function(file.kin0, file.kin=NULL, iids=NULL, type="kinship", verbose=TRUE){
	# read in the files
	if(verbose){ message("Reading .kin0 file...") }
	tmp0 <- read.table(file.kin0, header=TRUE);
	if(!is.null(file.kin)){
		if(verbose){ message("Reading .kin file...") }
		tmp1 <- read.table(file.kin, header=TRUE);
	}
	
	# find the unique IDs
	if(verbose){ message("Determining Unique Individual IDs from KING Output...") }
	uid01 <- as.character(unique(tmp0$ID1))
	uid02 <- as.character(unique(tmp0$ID2))
	if(!is.null(file.kin)){
		uid11 <- as.character(unique(tmp1$ID1))
		uid12 <- as.character(unique(tmp1$ID2))
		ids <- as.character(unique(c(uid01,uid02,uid11,uid12)))
	}else{
		ids <- as.character(unique(c(uid01,uid02)))
	}
	
	# check IDs
	if(!is.null(iids)){
		if(verbose){ message("Checking Provided Individual IDs") }
		if(!all(is.element(iids,ids))){
			stop("Some of the provided iids are not in the KING output")
		}
		if(!all(is.element(ids,iids))){
			stop("Some of the IDs in the KING output are not specified in the provided iids")
		}
		ids <- iids
	}
	
	# number of samples
	nsamp <- length(ids)
	
	# create empty matrix
	kingMat <- matrix(0, nrow=nsamp, ncol=nsamp)
	row.names(kingMat) <- ids; colnames(kingMat) <- ids;
	
	# add entries from .kin0
	if(verbose){ message("Adding Kinship Entries from .kin0 file...") }
	nr <- dim(tmp0)[1]
	id1 <- as.character(tmp0$ID1)
	id2 <- as.character(tmp0$ID2)
	for(r in 1:nr){
        if(type=="kinship"){
            kingMat[id1[r],id2[r]] <- tmp0$Kinship[r]
        }else if(type=="IBS0"){
            kingMat[id1[r],id2[r]] <- tmp0$IBS0[r]
        }
	}
	
	# add entries from .kin
	if(!is.null(file.kin)){
		if(verbose){ message("Adding Kinship Entries from .kin file...") }
		nr <- dim(tmp1)[1]
		id1 <- as.character(tmp1$ID1)
		id2 <- as.character(tmp1$ID2)
		for(r in 1:nr){
            if(type=="kinship"){
                kingMat[id1[r],id2[r]] <- tmp1$Kinship[r]
            }else if(type=="IBS0"){
                kingMat[id1[r],id2[r]] <- tmp1$IBS0[r]
            }
		}
	}
	
	# fill in lower half of the matrix
	kingMat <- kingMat + t(kingMat)
    diag(kingMat) <- 0.5
	
	# put ids as row and column names
	colnames(kingMat) <- ids
	rownames(kingMat) <- ids
	
	# return results
	return(kingMat)
	
}
