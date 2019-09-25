createClusterMatrix <- function(dat, id.colname, cluster.colname, ids.no.cluster = NULL, no.cluster.ids.pad = 0){
  stopifnot(is.element(no.cluster.ids.pad, c(0,1)))
  # in case these are factors: 
  dat$id <- dat[[id.colname]]
  dat$cluster.id <- dat[[cluster.colname]]
  
  dat <- dat[, c("id", "cluster.id")]
  dat[["id"]] <- as.character(dat[["id"]])
  dat[["cluster.id"]] <- as.character(dat[["cluster.id"]])
  
  if (!is.null(ids.no.cluster)){
    dat <- rbind(dat, data.frame(id = ids.no.cluster, cluster.id = ids.no.cluster, stringsAsFactors = FALSE))
  }
  
  cluster.matrix <- Diagonal(nrow(dat))
  colnames(cluster.matrix) <- rownames(cluster.matrix) <- dat[["id"]]
  
  if (!is.null(ids.no.cluster) & no.cluster.ids.pad != 1){
    diag(cluster.matrix)[match(ids.no.cluster, colnames(cluster.matrix))] <- no.cluster.ids.pad
  }
  
  
  # find all clusters with more than a single person:
  cluster.counts <- table(dat[["cluster.id"]])
  cluster.more.one <- names(cluster.counts[which(cluster.counts > 1)])
  
  for (i in seq_along(cluster.more.one)){
    cur.cluster <- cluster.more.one[i]
    person.ids.in.cluster <- dat[["id"]][which(dat[["cluster.id"]] == cur.cluster)]
    
    # make relevant entries in cluster.matrix equal to 1:
    cluster.matrix[person.ids.in.cluster, person.ids.in.cluster] <- 1
  }
  return(cluster.matrix)
}

