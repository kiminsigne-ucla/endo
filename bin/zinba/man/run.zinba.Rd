    \name{run.zinba}
          \Rdversion{1.0}
          \alias{run.zinba}
          \title{Run Analysis pipeline}
          \description{
This is the main function in Zinba consisting of three sequential steps:  
	(1) Processing mapped sample reads and calculating the correspoding covariate information, 
	(2) determining regions significantly enriched for reads given a set of covariates, and 
	(3) merging adjacent and overlapping windows and running a peak boundary refinement to get the exact boundaries of peak regions within these merged regions.  

The first step, the buildwindow step, is optional as one can elect to start their analysis using existing set of processed data by referencing its corresponding filelist file, defined below.  Also, the peak refinement step is also optional, where one can elect just to use the merged significant window coordinates instead of the refined boundaries (especially in the case of broader signal).  

If specified, the first step is taking the raw mapped sample reads, raw mapped input reads (if available), alignability directory, and current build of the genome and building the datasets needed to run the analysis to detect enriched regions.  The corresponding individual function to this is the buildwindowdata() function. After the data is built, the locations of each built file (one for each chromsome and offset) is placed in a filelist file.  

The next step uses the locations in filelist to import this built data into the analysis step to find windows that are likely to be enriched for counts given a set of covariates.  The corresponding individual function to this is getsigwindows(). The posterior probabilities for each window are saved, and the locations of these files are placed in the outfile.winlist file, where outfile is the name chosen to denote the output files of this current run.  The winlist is then fed to getrefinedpeaks() if peak refinement is desired. The overlapping signficant windows above the peakconfidence threshold are merged and the SBPC (basecount) for these merged windows are imported and exact peak boundaries are determined and outputted to outfile.peaks.   
          }
          \usage{run.zinba(buildwin=1,refinepeaks=1,outfile=NULL,peakconfidence=.99,twoBit=NULL,numProc=1,
    seq=NULL,input="none",filetype="bowtie",align=NULL,extension=NULL,selectmodel=NULL,
    #optional parameters below
    method='mixture',filelist=NULL,winSize=250,offset=0,cnvWinSize=100000, cnvOffset=0,
    formula=NULL,formulaE=NULL,formulaZ=NULL,threshold=0.01,basecountfile=NULL,tol=10^-5,
    winGap=0,pWinSize=200,pquant=1,printFullOut=0,initmethod='count',diff=0,
    cleanup=FALSE,selectchr=NULL, selecttype="dirty", selectcovs=NULL

)
    



	}					
          \arguments{
	    \item{buildwin}{Default value is 1, specifying to build the window files from raw data sources.  Raw data corresponds to the mapped sample reads, raw mapped input reads (if available), the current build of the human genome and the alignabilty directory.  A value of 0 specifies you would like to skip this step and use an existing filelist (corresponding to a previously built set of window data)}
            \item{refinepeaks}{Whether exact boundaries of peaks within merged significantly enriched regions is requested (default), otherwise  set to 0.  Helps to isolate punctuate, sharp peaks.  If TRUE, then specify parameter basecountfile}
	    \item{outfile}{Prefix used to denote the .wins, .winlist and .peaks files that are outputted by zinba}
	    \item{peakconfidence}{Posterior probability threshold for signfificant windows if method='mixture'.  Higher is more stringent, must be between 0 and 1. Must be specified if method=mixture, if method='pscl' then specify threshold instead.}
	    \item{twoBit}{Path to build of the genome your reads were mapped to, in .2bit format}
	    \item{numProc}{Number of concurrent jobs to run in parallel, default is 1, recommended to use more if more computing cores are available.  More memory is needed however when more processes are run in parallel.  If using winSize of 500, 4 cores is sufficient for a desktop with 4 GB of memory, smaller winSizes require more memory}
	    \item{seq}{Path to mapped sample reads if buildwin=1, formatted as either 'bed', 'tagAlign', or 'bowtie'}
	    \item{input}{Path to mapped input reads if buildwin=1, formatted as either 'bed', 'tagAlign', or 'bowtie'.  If left blank, then defaults to 'none'.  Input control is not necessary to run ZINBA}
	    \item{filetype}{Format of mapped sample and input reads.  'bed' is files in the standard .bed format, 'tagAlign' signifies those in .taf format, and 'bowtie' signifies mapped reads directly outputted from bowtie.  Default is 'bowtie'}
	    \item{align}{Path to directory containing alignability files for each chromsome, obtained from alignAdjust (check how these files are generated) or downloadd from our respository.  Alignability information is specific to the uniquness threshold one used to initially filter their mapped reads and the length of the sequence tags used}
	    \item{extension}{Average length of fragments in fragment library used, typically around 200}
	    \item{selectmodel}{if TRUE, then perform model selection before mixture regression to determine best set of component covariates.  If TRUE parameters selectchr,selecttype,and selectcovs need to be specified.  Although time consiming it is suggested to be done the first time you are running a new or unfamilar type of data, otherwise specify formula, formulaE, and formulaZ.}   
	    \item{method}{OPTIONAL. If 'mixture' is selected, the mixture regression model approach is used (must specify formula and formulaE). If 'pscl' is specified, then only needs to specify formula.  Pscl is an adhoc method that is used on FAIRE data}
	   \item{threshold}{FDR threshold for significant windows for pscl method, lower is more stringent.  Must be specified if method='pscl', otherwise ignore.}
	    \item{filelist}{OPTIONAL. If buildwin=0, this parameter is a path to an existing filelist corresponding to previously built dataset (so data does not have to be built again for analysis)}
	    \item{winSize}{OPTIONAL.  Window size for build window data, default=250 bp}
	    \item{offset}{OPTIONAL. Bp offset, default is no offsets (offset=0).  If one's winSize is 500 and would like 4 equally spaced offsets, then specifying offset=125 would achieve this}
	    \item{cnvWinSize}{OPTIONAL. Size of windows used to calculate CNV activity in sample, default size is 100000}
	    \item{cnvOffset}{OPTIONAL. Offset for CNV windows, typically 2500bp is suffcient, although default is no offsets}
	    \item{formula}{Specifies the formula desired to model the background component of counts.  Must begin with "exp_count~"  (without quotes) and end with a set of desired covariates separated by "+" or "*" (without quotes). The covariates could be input_count (input), gcPerc (gc percentage), exp_cnvwin_log (cnv estimate), and align_perc (alignability estimate).  The "+" indicates linear addition of terms, and "*" indicates linear addition of terms plus an interaction between the terms.  Examples:  If sufficiently sequenced input is present, exp_count~input_count.  If no input is present, exp_count~gcPerc+exp_cnvwin_log+align_perc is sufficient. }
	    \item{formulaE}{Only used if method is "mixture".  Specifies the formula desired to model the enriched component.  If one is unsure about making assumptions about the relationship with covarites and enriched signal, use exp_count~1}
	    \item{formulaZ}{Only used if method is "mixture".  Specifies the formula desired to model the zero inflated component.  If one is unsure about making assumptions about the relationship with covarites and enriched signal, leave blank and will use the model specified for background (formula)}	
	    \item{selectchr}{chromsome to use for model selection typicaly a smaller one such as human chromosome 22 (selectchr="chr22")}
	    \item{selectcovs}{Set of starting covariates to consider for model selection entered as a vector of strings representing each covariate to be used.  For example, to only consider covariates GC content, alignability, and input control then enter c("gcPerc", "align_perc", "input_count").  To only consider covariates GC content, alignability, and local background then enter c("gcPerc", "align_perc", "exp_cnvwin_log").  It is advised not to include the local background estimate in model selection unless input is not available.}
	    \item{selecttype}{OPTIONAL.  If "dirty" (default) then covariates for background and enriched components are determined while zero-inflation formula is held fixed with no covariates.  After background and enriched covariates are determined, the best zero inflation formula is determined given the best background and enriched covariates. If set to "complete" then all possible sets of covariates are determined (guaruntees best model will be selected but 20x longer computation time).  This option should be used if estimates of covariate effects are desired rather than just peak calling}
	    \item{initmethod}{OPTIONAL.  Used only if method is "mixture".  Method of intial partitioning of windows into each component prior to mixture regression model run.  usually does nto make a difference. 'count' is default, specifies taking the windows with largest counts and classifying them intially as enriched.  'quantile' runs a quantile regression onthe data using the background covariates, takes the top residuals and assigns them as enriched.  'pscl' uses the pscl method to find likely enriched regions and assigns these windows to the enriched component. } 
	    \item{printFullOut}{OPTIONAL. If set to 1, prints out the original dataset along with the posterior probabilities.  Otherwise, prints out only window coordinates and significance score (default)}
	    \item{tol}{OPTIONAL. Relative tolerance for convergence, default  is 10^-5}
	    \item{basecountfile}{OPTIONAL.  Must be specified if refinepeaks is 1.  Path to basecount track containing SBPC information for the entire genome, generated by basealigncount}
	    \item{pWinSize}{OPTIONAL. In the peak refinement step, reduce this parameter to get more sensitivity to detect local Maximums, default is 200.  This may increase false positives}
	    \item{pquant}{OPTIONAL. In the peak refinement step, filters out local maximums whose height is lower than the 'pquant' quantile in the window. Use pquant=1 for selecting only the global max in a collapsed region, god for ChIP-seq datasets.  Default is .75}
	   \item{diff}{OPTIONAL.  Unimplemented parameter, pertains to differential peak expression analysis accross samples}	
   	   \item{cleanup}{OPTIONAL. If set equal to TRUE, will delete intermediate output file folder.  Default value is FALSE.}
}
          \seealso{
            \code{\link{save}}.
          }
          \examples{
#Chip-Seq with recommended parameters
run.zinba(
   buildwin=1,
   refinepeaks=1,
   method='mixture', 
   peakconfidence=.99,
   outfile='data/ctcf',
   twoBit='hg18.2bit',
   numProc=4,
   seq='data/ctcfGM12878rep3chr22.taf' 
   filelist='data/ctcf.list',
   input='data/inputGM12878rep3chr22.taf'
   filetype='tagAlign', 
   align='align1/',
   winSize=500,
   offset=125,
   cnvWinSize=100000,
   cnvOffset=2500,
   extension=200,
   formula=exp_count~input_count,
   formulaE=exp_count~1,
   formulaZ=exp_count~1,
   initmethod='count', 
   printFullOut=1,
   tol=10^-5,
   pWinSize=200,
   pquant=1,
   basecountfile='ctcfGM12878rep3chr22.basecount',
   diff=0,
   cleanup=FALSE
   )	   
		
#FAIRE-seq with recommended parameters, no input.  
run.zinba(
   buildwin=1,
   refinepeaks=1,
   method='mixture', 
   peakconfidence=.95,
   outfile='data/faire',
   twoBit='hg18.2bit',
   numProc=4,
   seq='data/faireGM12878rep1chr22.taf',
   filelist='data/faire.list',
   input='none',
   filetype='tagAlign', 
   align='align4/',
   winSize=250,
   offset=50,
   cnvWinSize=100000,
   cnvOffset=2500,
   extension=134,
   formula=exp_count~exp_cnvwin_log+gcPerc+align_perc,
   formulaE=exp_count~1,
   formulaZ=exp_count~exp_cnvwin_log+gcPerc+align_perc,
   initmethod='count', 
   printFullOut=1,
   tol=10^-5,
   pWinSize=200,
   pquant=1,
   basecountfile='data/faireGM12878rep1chr22.basecount',
   diff=0,
   cleanup=FALSE
   )



   
          }
          \keyword{file}
