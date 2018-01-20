    \name{alignAdjust}
          \Rdversion{1.0}
          \alias{alignAdjust}
          \title{Custom Alignability Track Generation}
          \description{
            Part of generateAlignability.  From the peakseq mappability file for a given sequence read length and genome, alignAdjust assigns a score of 0, 1, or 2 to every base pair location in the genome.  A score of 2 indicates that the length-n sequence read starting at that particular base meets the alignability threshold specified in either strand direction, 1 if it only meets this threshold for one strand direction, and 0 if the the length-n sequence in either strand direction meets this threshold. To account for the fact that reads are extended to estimated length of the fragments used in the sequence library, we shift the scores by this adjustment length as well so that a reads counted in a given window will have the same matching alignability criteria. 
          }
          \usage{alignAdjust(inputfile, outDir="",twoBitFile,athresh=1,adjustsize=0)}
          \arguments{
	    \item{inputfile}{file containing a score at every base corresponding to the number of times a given k-mer occurs, originating from that base; data is formatted as a fixedStep wiggle track ( see
UCSC Genome Browser custom tracks).  This file is outputted from convertmappability}
	    \item{outDir}{directory where results will be stored}
	    \item{twobitfile}{path to current build of human genome in .2bit format}
	    \item{athresh}{uniqueness threshold, number of occurrences of a given k-mer imposed during alignment (1 means absolute uniqueness)}
	    \item{adjustsize}{adjustsize one half of the extension length used for aligned reads}
	  }
          \seealso{
            \code{\link{save}}.
          }
          \examples{
#Corresponding to FAIRE example in run.zinba.  For FAIRE, four or less matches to the 
    genome allowed, so athresh=4 rather than 1
alignAdjust('mappability.wig', outDir="align/",twoBitFile,athresh=4,adjustsize=134/2)
   
          }
          \keyword{file}
