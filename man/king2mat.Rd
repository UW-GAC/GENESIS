\name{king2mat}
\alias{king2mat}
\title{Convert KING text output to an R Matrix}
\description{\code{king2mat} is used to extract the pairwise kinship coefficient estimates or IBS0 values from the output text files of KING and put them into an R object of class \code{matrix} that can be read by the functions \code{\link{pcair}} and \code{\link{pcairPartition}}.}
\usage{
king2mat(file.kin0, file.kin = NULL, iids = NULL, 
            type = "kinship", verbose = TRUE)
}
\arguments{
  \item{file.kin0}{File name of the .kin0 text file output from KING.}
  \item{file.kin}{Optional file name of the .kin text file output from KING.}
  \item{iids}{An optional vector of individual IDs in the same order as desired for the output matrix.  See 'Details' for more information.}
  \item{type}{Character string taking the values "kinship" (default) or "IBS0", to inform the function to read in kinship coefficeints or IBS0 values from the KING output.}
  \item{verbose}{A logical indicating whether or not to print status updates to the console; the default is TRUE.}
}
\details{When using the function \code{pcair}, it is important that the order of individuals in the \code{kinMat} matrix matches the order of individuals in \code{genoData}.  The KING software has a tendency to reorder individuals.  If \code{iids = NULL}, the default is for the order to be taken from the KING output text file.  By specifying \code{iids} the user can control the order of individuals in the output matrix.  The IDs used for \code{iids} must be the same set of character IDs that are output as columns 'ID1' and 'ID2' in the KING output text files; all of the IDs specified in \code{iids} must be in the KING output, and all IDs in the KING output must be specified in \code{iids}.}
\value{An object of class '\code{matrix}' with pairwise kinship coefficients or IBS0 values as estimated by KING for each pair of individuals in the sample.  The estimates are on both the upper and lower triangle of the matrix, and the diagonal is arbitrailly set to 0.5.  Individual IDs are set as the column and row names of the matrix.}
\references{
    Conomos M.P., Miller M., & Thornton T. (2015). Robust Inference of Population Structure for Ancestry Prediction and Correction of Stratification in the Presence of Relatedness. Genetic Epidemiology, 39(4), 276-293.
    
    Manichaikul, A., Mychaleckyj, J.C., Rich, S.S., Daly, K., Sale, M., & Chen, W.M. (2010). Robust relationship inference in genome-wide association studies. Bioinformatics, 26(22), 2867-2873.
}
\author{Matthew P. Conomos}
\seealso{
	\code{\link{pcair}} and \code{\link{pcairPartition}} for functions that use the output matrix.
}
\examples{
file.kin0 <- system.file("extdata", "MXL_ASW.kin0", package="GENESIS")
file.kin <- system.file("extdata", "MXL_ASW.kin", package="GENESIS")
KINGmat <- king2mat(file.kin0 = file.kin0, file.kin = file.kin, type="kinship")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
