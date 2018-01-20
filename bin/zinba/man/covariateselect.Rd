    \name{covariateselect}
          \Rdversion{1.0}
          \alias{covariateselect}
          \title{Determine Best Sets of Covariates to Model Each Component}
          \description{
            Implements model selection procedure.  Getsigwindows() is called using modelselect=TRUE to determine the best fitting model on a subset of the data.  Each model's BIC is calculated and is outputted to a file at location "loc". The model with the lowest score is then chosen.    All possible models given a set of starting covariates are considered when selection="complete".  A much faster model selection procedure (selection="dirty") fixes the zero-inflation component to have no covariates and first selects best covariates for the background and enriched components.  After these covariates are selected, the best zero-inflation covariates are chosen.  While this procedure is quicker, it is not guarunteed that the best model will be chosen and should not be used if one cares more about selecting the best model and estimated covariate effects rather than peak calling alone. 
          }
          \usage{covariateselect(file=NULL, covs=NULL, selection="dirty",formulaZ=NULL,numProc=1,loc=NULL, start=1, append=FALSE)}
          \arguments{
	    \item{file}{Path to an existing built data file corresponding to a particular chromsome and offset, typically smaller chromosome}
	    \item{selection}{If "dirty" (default) then covariates for background and enriched components are determined while zero-inflation formula is held fixed with no covariates.  After background and enriched covariates are determined, the best zero inflation formula is determined given the best background and enriched covariates. If set to "complete" then all possible sets of covariates are determined (guaruntees best model will be selected but 20x longer computation time).  This option should be used if estimates of covariate effects are desired rather than just peak calling}
	    \item{covs}{Set of starting covariates to consider for model selection entered as a vector of strings representing each covariate to be used.  For example, to only consider covariates GC content, alignability, and input control then enter c("gcPerc", "align_perc", "input_count").  To only consider covariates GC content, alignability, and local background then enter c("gcPerc", "align_perc", "exp_cnvwin_log").  It is advised not to include the local background estimate in model selection unless input is not available.}
	    \item{formulaZ}{OPTIONAL.  Formula to fix the zero-inflation component to if selection="dirty".  Default is no covariates}
	    \item{numProc}{OPTIONAL.  Number of processors to use}
	    \item{start}{OPTIONAL.  Start the model at a certain model number rather from just the beginning}
	    \item{append}{OPTIONAL. append model selection results to an existing model file}
          }
	  \seealso{
            \code{\link{save}}.
          }
          \examples{
#Corresponding to FAIRE example in run.zinba:
covariateselect(file="FAIRE/faire_chr10_win250bp_offset0bp.txt;FAIRE/faire_chr10_win250bp_offset50bp.txt", 
		covs=c("gcPerc","align_perc"), 
		selection="dirty",
		numProc=1,
		loc="data.model")
          }
          \keyword{file}
