\name{sample_annotation_1KG}
\alias{sample_annotation_1KG}
\docType{data}
\title{Annotation for 1000 genomes Phase 3 samples}
\description{
    Annotation for 1000 genomes Phase 3 samples included in the VCF files in "extdata/1KG".
}
\usage{
    data(sample_annotation_1KG)
}
\format{
  A data.frame with columns:
  \itemize{
    \item{sample.id}{Sample identifier}
    \item{Population}{Population of sample}
    \item{sex}{Sex of sample}
  }
}
\source{
	ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp
}
\references{
	A global reference for human genetic variation, The 1000 Genomes Project Consortium, Nature 526, 68-74 (01 October 2015) doi:10.1038/nature15393.
}
\keyword{datasets}
