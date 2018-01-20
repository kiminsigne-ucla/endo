    \name{covariateselect}
          \Rdversion{1.0}
          \alias{covariateselect}
          \title{Determine Best Sets of Covariates to Model Each Component}
          \description{
            Implements model selection procedure.  Getsigwindows() is called using modelselect=TRUE to determine the best fitting model on a subset of the data.  Each model's BIC is calculated and is outputted to a file at location "loc". The model with the lowest score is then chosen.    All possible models given a set of starting covariates are considered when selecttype="complete".  A much faster model selection procedure (selecttype="dirty") fixes the zero-inflation component to have no covariates and first selects best covariates for the background and enriched components.  After these covariates are selected, the best zero-inflation covariates are chosen. 

          }
          \usage{covariateselect(inputfile, outputfile)}
          \arguments{
	    \item{inputfile}{Path to file in the current directory map.list, which contains the names of the chromosome specific binary files in the mappability directory. Each chromosome binary file contains a score at every base corresponding to the number of times a given k-mer occurs in the genome, originating from that base; downloaded from repository or self generated}
	    \item{outputfile}{single wig file where results should be printed}
          }
	  \seealso{
            \code{\link{save}}.
          }
          \examples{
#Corresponding to FAIRE example in run.zinba:
covariateselect(inputfile='mappabilitylist.txt',  outputfile='mappability.wig')
   
          }
          \keyword{file}
