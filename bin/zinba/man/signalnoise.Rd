    \name{signalnoise}
          \Rdversion{1.0}
          \alias{twobittofa}
          \title{Convert Fasta Formatted Files From .2bit}
          \description{
            Computes Signal to Noise information for a given basecountfile.  The log ratio of the max to median basecount in every window is computed, and the summary information including the mean and median log ratio is given.  Higher mean values of this log ratio is a crude measure of better signal to noise.  
          }
          \usage{signalnoise(inputfile,twoBitFile,winSize=100000)}
          \arguments{
	    \item{inputfile}{Path to basecountfile to be analyzed}
	    \item{twoBitFile}{Path to corresponding genome in .2bit format}
	    \item{winSize}{Window Size for analysis, default is 100kb}
	  }
          \seealso{
            \code{\link{save}}.
          }
          \examples{
		signalnoise('data/faireGM12878rep1chr22.basecount', 'hg18.2bit', winSize=5000)
		
          }
          \keyword{file}
