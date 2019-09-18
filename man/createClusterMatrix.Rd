\name{createClusterMatrix}
\alias{createClusterMatrix }

\title{Create matrix representation of the associations of a set of individuals with clusters}
\description{
Creates a sparse matrix tabulating the cluster associations of individuals in a data.frame, and potentially also an additional set of individuals that are not associated with any cluster. The i,j entry in the matrix will have the value 1 if individuals i and j are in the same cluster, and 0 otherwise.
}
\usage{
createClusterMatrix(dat, id.colname, cluster.colname, ids.no.cluster = NULL, no.cluster.ids.pad = 0)
}
\arguments{
    \item{dat}{A data.frame containing a column with individual/sample (different for each matrix row) IDs, and column with cluster IDs.}
    \item{id.colname}{Name of the column with individual/sample IDs}
    \item{cluster.colname}{The name of the column with cluster IDs}
    \item{ids.no.cluster}{An optional parameter. A vector of additional IDs not in the \code{data.frame} dat, that are not associated with the clusters provided for the individuals in dat}
    \item{no.cluster.ids.pad}{An optional parameter. If there are IDs provided in \code{ids.no.cluster}, such that they do not have any cluster specified, matrix tabulating cluster associations will have no.cluster.ids.path number on its diagonal elements corresponding to these individuals.  Default number is \code{0}}
}
\author{Tamar Sofer}

\examples{
	n <- 200
	n.cluster <- 20
	cluster.names <- paste0("cluster.", 1:n.cluster)
	# in this toy data set, the first 100 people are from 20 clusters. The last 100 people do not belong to any pre-defined cluster, and therefore are defined to have their own, one-person, cluster.
	dat <- data.frame(person.id = paste0("person.", 1:n), 
	                  cluster.id = c(sample(cluster.names, n/2, replace = TRUE), paste0("person.", (n/2 + 1):n)),
	                  stringsAsFactors = FALSE)
	 
	cluster.mat <- createClusterMatrix(dat, "person.id", "cluster.id")                                             
	 
	 
	## now assume we're comgining two studies (e.g. HCHS/SOL with household IDs, and WHI):
	 
	dat <- data.frame(person.id = paste0("person.", 1:n), 
	                  cluster.id = sample(cluster.names, n, replace = TRUE),
	                  stringsAsFactors = FALSE)
	 
	ids.no.cluster  <-              paste0("person.s2.", 1:n)                                           
	cluster.mat <- createClusterMatrix(dat, "person.id", "cluster.id", ids.no.cluster = ids.no.cluster, no.cluster.ids.pad = 0) 
}