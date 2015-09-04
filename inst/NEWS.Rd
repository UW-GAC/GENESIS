\name{NEWS}
\title{NEWS for GENESIS}
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
