    \name{getsigwindows}
          \Rdversion{1.0}
          \alias{getsigwindows}
          \title{Determine Significantly Enriched Windows}
          \description{
            Determines the set of enriched windows using either a residual-based or mixture model
          }
          \usage{getsigwindows(file, formula, formulaE, formulaZ, threshold=.01, peakconfidence=.99, 
    winout, printFullOut=0, tol=10^-5, method='mixture', initmethod="count", 
    diff=0)}
          \arguments{
	    \item{file}{Path to a build window data file for a chomosome containing window counts and covariate information.  Specifying a string of offsets for a chromosome separated by a ";" will apply get sigwindows for each path in the string and output it all results in one file at the specified winout path.  In run.zinba, these paths are read from a .list filelist file where offsets are separated by  “;”.  }
	    \item{winout}{Path to where results from mixture regression are outputted.  The window ID (coordinates pastable into the UCSC Genome Browser, original data matrix, filtering status, and posterior enrichment probability are printed}
            \item{method}{if 'mixture' is selected, the mixture regression model approach is used . If 'pscl' is specified, then only needs to specify formula.  Pscl is an adhoc method that is used on FAIRE data, not actively maintained}
	    \item{formula}{Specifies the formula desired to model the background component of counts.  Must begin with "exp_count~"  (without quotes) and end with a set of desired covariates separated by "+" or "*" (without quotes). The covariates could be input_count (input), gcPerc (gc percentage), exp_cnvwin_log (cnv estimate), and align_perc (alignability estimate).  The "+" indicates linear addition of terms, and "*" indicates linear addition of terms plus an interaction between the terms.  Examples:  If sufficiently sequenced input is present, exp_count~input_count.  If no input is present, exp_count~gcPerc+exp_cnvwin_log+align_perc is sufficient. }
	    \item{formulaE}{Only used if method is "mixture".  Specifies the formula desired to model the enriched component.  If one is unsure about making assumptions about the relationship with covarites and enriched signal, use exp_count~1}
	    \item{formulaZ}{Only used if method is "mixture".  Specifies the formula desired to model the zero inflated component.  If one is unsure about making assumptions about the relationship with covarites and enriched signal, leave blank and will use the model specified for background (formula)}
	    \item{peakconfidence}{posterior probability threshold  for signfiicant windows if the mixture regression model is used}
	    \item{threshold}{FDR threshold for signfiicant windows if the pscl method is used}
	    \item{tol}{OPTIONAL.  Relative tolerance for convergence, default  is 10^-5}
	    \item{initmethod}{OPTIONAL. Method of intial partitioning of windows into each component prior to mixture regression model run.  usually does nto make a difference. 'count' is default, specifies taking the windows with largest counts and classifying them intially as enriched.  'quantile' runs a quantile regression onthe data using the background covariates, takes the top residuals and assigns them as enriched.  'pscl' uses the pscl method to find likely enriched regions and assigns these windows to the enriched component. } 
	    \item{printFullOut}{OPTIONAL.  If set to 1, prints out the original dataset along with the posterior probabilities (default).  Otherwise, prints out only window coordinates and posterior probabilites (uses less disk space but  cannot explore the data's basic relationships with posterior probability and covariates later if desired}
	 \item{modelselect}{OPTIONAL. Internal parameter to notify getsigwindows that it is being used for model selection, not peak calling.  Users will not need to typically use this}
}
          \seealso{
            \code{\link{save}}.
          }
          \examples{
	   #FAIRE-seq with recommended parameters, no input, multiple offsets.  
run.zinba(
    file='FAIRE/faire_chr10_win250bp_offset0bp.txt;FAIRE/faire_chr10_win250bp_offset50bp.txt', 
    winout='FAIRE/faire', formula=exp_count~exp_cnvwin_log+gcPerc+align_perc, formulaE=exp_count~1,  formulaZ=exp_count~exp_cnvwin_log+gcPerc+align_perc, peakconfidence=.99, method='mixture', )


   
          }
          \keyword{file}
