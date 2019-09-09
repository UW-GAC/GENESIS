\name{makeSparseMatrix}
\alias{makeSparseMatrix}
\alias{makeSparseMatrix-methods}
\alias{makeSparseMatrix,matrix-method}
\alias{makeSparseMatrix,Matrix-method}
\alias{makeSparseMatrix,data.frame-method}
\alias{makeSparseMatrix,data.table-method}
\title{Make a sparse matrix from a dense matrix or a table of pairwise values}
\description{\code{makeSparseMatrix} is used to create a sparse block-diagonal matrix from a dense matrix or a table of pairwise values, using a user-specified threshold for sparsity. An object of class \code{Matrix} will be returned. A sparse matrix may be useful for fitting the association test null model with \code{\link{fitNullModel}} when working with very large sample sizes.}
\usage{
\S4method{makeSparseMatrix}{data.table}(x, thresh = NULL, sample.include = NULL, diag.value = NULL,
                 verbose = TRUE)
\S4method{makeSparseMatrix}{data.frame}(x, thresh = NULL, sample.include = NULL, diag.value = NULL, 
                 verbose = TRUE)
\S4method{makeSparseMatrix}{matrix}(x, thresh = 2^(-11/2), sample.include = NULL, diag.value = NULL, 
                 verbose = TRUE)
\S4method{makeSparseMatrix}{Matrix}(x, thresh = 2^(-11/2), sample.include = NULL, diag.value = NULL, 
                 verbose = TRUE)
}
\arguments{
  \item{x}{An object to coerce to a sparse matrix. May be of class \code{matrix}, \code{Matrix}, \code{data.frame}, or \code{data.table}. When a \code{matrix} or \code{Matrix}, row and column names should be set to sample.ids; when a data.frame or data.table, should have 3 columns: ID1, ID2, and value.}
  \item{thresh}{Value threshold for clustering samples to make the output matrix sparse block-diagonal. When \code{NULL}, no clustering is done. See 'Details'.}
  \item{sample.include}{An optional vector of sample.id indicating all samples that should be included in the output matrix; see 'Details' for usage.}
  \item{diag.value}{When \code{NULL} (by Default), values for the diagonal of the output matrix should be provided in \code{x}. Setting \code{diag.value} to a numeric value will make all values on the diagonal of the output matrix that value.}
  \item{verbose}{A logical indicating whether or not to print status updates to the console; the default is TRUE.}
}
\details{\code{sample.include} has two primary functions: 1) It can be used to subset samples provided in \code{x}. 2) \code{sample.include} can include sample.id not in \code{x}; these additional samples will be included in the output matrix, with a value of 0 for all off-diagonal elements, and the value provided by \code{diag.value} for the diagonal elements. When left \code{NULL}, the function will determine the list of samples from what is observed in \code{x}.

\code{thresh} sets a threhsold for clustering samples such that any pair with a value greater than \code{thresh} is in the same cluster. All values within a cluster are kept, even if they are below \code{thresh}. All values between clusters are set to 0, creating a sparse, block-diagonal matrix. When \code{thresh} is \code{NULL}, no clustering is done and all samples are returned in one block. This feature is useful when converting a data.frame of kinship estimates to a matrix.}
\value{An object of class '\code{Matrix}'. Samples may be in a different order than in the input \code{x}, as samples are sorted by ID or rowname within each block (including within the block of unrelateds).}

\author{Matthew P. Conomos}
\seealso{
	\code{\link{kingToMatrix}} and \code{\link{pcrelateToMatrix}} for functions that use this function.
}
% \examples{
% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
