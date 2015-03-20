\name{plot.pcair}
\alias{plot.pcair}
\title{PC-AiR: Plotting PCs}
\description{\code{plot.pcair} is used to plot pairs of principal components contained in a class '\code{pcair}' object obtained as output from the \code{pcair} function.}
\usage{
\method{plot}{pcair}(x, vx = 1, vy = 2, pch = NULL, col = NULL, 
        xlim = NULL, ylim = NULL, main = NULL, sub = NULL, 
        xlab = NULL, ylab = NULL, ...)
}
\arguments{
  \item{x}{An object of class '\code{pcair}' obtained as output from the \code{pcair} function.}
  \item{vx}{An integer indicating which principal component to plot on the x-axis; the default is 1.}
  \item{vy}{An integer indicating which principal component to plot on the y-axis; the default is 2.}
  \item{pch}{Either an integer specifying a symbol or a single character to be used in plotting points.  If \code{NULL}, the default is dots for the 'unrelated subset' and + for the 'related subset'.}
  \item{col}{A specification for the plotting color for points.  If \code{NULL}, the default is black for the 'unrelated subset' and blue for the 'related subset'.}
  \item{xlim}{The range of values shown on the x-axis.  If \code{NULL}, the default shows all points.}
  \item{ylim}{The range of values shown on the y-axis.  If \code{NULL}, the default shows all points.}
  \item{main}{An overall title for the plot.  If \code{NULL}, the default specifies which PC-AiR PCs are plotted.}
  \item{sub}{A sub title for the plot.  If \code{NULL}, the default is none.}
  \item{xlab}{A title for the x-axis.  If \code{NULL}, the default specifies which PC-AiR PC is plotted.}
  \item{ylab}{A title for the y-axis.  If \code{NULL}, the default specifies which PC-AiR PC is plotted.}
  \item{\dots}{Other parameters to be passsed through to plotting functions, (see \code{\link{par}}).}
}
\details{
	This function provides a quick and easy way to plot principal components obtained with the function \code{pcair} to visualize the population structure captured by PC-AiR.
}
\value{A figure showing the selected principal components plotted against each other.}
\author{Matthew P. Conomos}
\seealso{
	\code{\link{pcair}} for obtaining principal components that capture population structure in the presence of relatedness.
	\code{\link{par}} for more in depth descriptions of plotting parameters.
	The generic function \code{\link{plot}}.
}
\examples{
# file path to GDS file
gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
# read in GDS data
HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
# create a GenotypeData class object
HapMap_genoData <- GenotypeData(HapMap_geno)
# load saved matrix of KING-robust estimates
data("HapMap_ASW_MXL_KINGmat")
# run PC-AiR
mypcair <- pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, 
                divMat = HapMap_ASW_MXL_KINGmat)
# plot top 2 PCs
plot(mypcair)
# plot PCs 3 and 4
plot(mypcair, vx = 3, vy = 4)
close(HapMap_genoData)
}

