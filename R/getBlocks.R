getBlocks <- function(n, block.size){
    # number of blocks
    nblocks <- ceiling(n/block.size)
    # start and end positions for blocks
    if(nblocks==1){
        block.start <- 1
        block.end <- n
    }else{
        block.start <- (0:(nblocks-1))*block.size+1
        block.end <- c( (1:(nblocks-1))*block.size, n )
    }

    return(list(start = block.start, end = block.end, n=nblocks))
}