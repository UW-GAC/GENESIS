\name{NEWS}
\title{NEWS for GENESIS}
\section{Version 2.11.4}{
    \itemize{
        \item Added assocTestSingle and fitNullModel methods for GenotypeData objects.
    }
}
\section{Version 2.9.3}{
    \itemize{
	\item New methods assocTestSingle and assocTestAggregate are refactors of assocTestMM and assocTestSeq/assocTestSeqWindow, respectively. assocTestSeq and assocTestSeqWindow are deprecated. assocTestMM is still used for GenotypeData objects, but will be deprecated in a future release. fitNullModel is a refactor of fitNullMM/fitNullReg and should be used with the new association test methods.
    }
}
\section{Version 2.7.4}{
    \itemize{
	\item In fitNullMM, for Binomial and Poisson GLM families, the variance of the fixed effects will not be multiplied by the dispersion parameter RSS.
    }
}
\section{Version 2.7.3}{
    \itemize{
	\item Change defaults in assocTestSeq*: Davies method for SKAT p-values, flat weights.
    }
}
\section{Version 2.6.0}{
    \itemize{
	\item Major bug fix: assocTestSeq no longer drops some variants from aggregate tests in the case where the same variants are included in more than one aggregate unit.
        \item Added function for analysis of admixture mapping data.
    }
}
\section{Version 2.2.0}{
    \itemize{
	\item Added support for sequence data files created with the SeqArray package.
        \item Added functions for Burden and SKAT tests of sequence data.
    }
}
\section{Version 2.0.0}{
    \itemize{
        \item Added functions for PC-Relate.  PC-Relate provides model-free estimation of recent genetic relatedness in general samples.  It can be used to estimate kinship coefficients, IBD sharing probabilities, and inbreeding coefficients using genome-wide SNP data.  PC-Relate accounts for population structure (ancestry) among sample individuals through the use of ancestry representative principal components (PCs) to provide accurate relatedness estimates due only to recent family (pedigree) structure.
        \item GENESIS now imports the package gdsfmt.
    }
}
\section{Version 1.1.1}{
	\itemize{
		\item Changed a few default values:  MAF filter from 0.05 to 0.01; kin.thresh from 0.025 to 2^(-11/2) ~ 0.022; div.thresh from -0.025 to -2^(-11/2) ~ -0.022; v (the number of PCs returned) from 10 to 20
		\item Changed the argument 'block.size' to 'snp.block.size'
		\item Slightly altered the command line reporting when verbose = TRUE
	}
}
\section{Version 0.99.4}{
    \itemize{
        \item Fixed a minor bug to zero out diagonal of divMat in pcairPartition()
    }
}
\section{Version 0.99.0}{
    \itemize{
        \item Initial version of GENESIS.  Contains functionality to perform PC-AiR.
    }
}
