    \name{basealigncount}
          \Rdversion{1.0}
          \alias{basealigncount}
          \title{Get SBC (basecount) for Mapped Sample Reads}
          \description{
            This function creates the SBC (basecount) file required by run.zinba that used for boundary detection in the peak refinement step.  This should be generated prior to running run.zinba if refinement is desired.   
          }
          \usage{basealigncount(inputfile,outputfile,twoBitFile,extension,filetype="bowtie")}
          \arguments{
	    \item{inputfile}{path to mapped sample reads, formatted as either 'bed', 'tagAlign', or 'bowtie'}
	    \item{outputfile}{path to file where basecount file will be outputted}
	    \item{twobitfile}{path to current build of the genome your reads where mapped to, in .2bit format}
	    \item{extension}{average length of fragments in fragment library used, typically 150-300 bp}
	    \item{filetype}{format of mapped sample and input reads.  'bed' is files in the standard .bed format, 'tagAlign' signifies those in .taf format, and 'bowtie' signifies mapped reads directly outputted from bowtie.  Default is 'bowtie'}
	  }
          \seealso{
            \code{\link{save}}.
          }
          \examples{
	   #Corresponding to FAIRE example in run.zinba:

basealigncount(
   inputfile='data/faireGM12878rep1chr22.taf', 
   outputfile='data/faireGM12878rep1chr22.basecount', 
   extension=134, 
   filetype='tagAlign', 
   twoBitFile='hg18.2bit' )
   
          }
          \keyword{file}
