    \name{buildwindowdata}
          \Rdversion{1.0}
          \alias{buildwindowdata}
          \title{Build Datasets for Mixture Regression Model}
          \description{
            Builds the datasets needed to run ZINBA's mixture regression model from raw data sources.  Takes the raw mapped sample reads, raw mapped input reads (if available), alignability directory, and current build of the human genome to calculate GC proportion, proportion of alignable bases, matching input control read counts, and a local background estimate.  These values are calculated for non-overlapping windows of specified size spanning each chromsome.  If a window offset is specified, additional data files are generated for windows specified at the offset distance. For example, if 500 bp windows and an offset of 125 bp are specified, then 4 offset files are created at 0bp, 125bp, 250bp, and 375 bp offsets.  In run.zinba, these files are automatically placed in the outfile_files directory, where "outfile" is the outfile prefix specified.  
          }
          \usage{buildwindowdata(seq, input="none", align, twoBit, winSize=250, offset=0, 
    cnvWinSize=100000, cnvOffset=2500, filelist, filetype="bowtie", 
    extension, outdir)}
          \arguments{
            \item{filelist}{At path 'filelist', a file will be created holding the paths to each of the built datasets.  Each line within this file pertains to each chromsome, and within each line chromosome offsets are separated by ";".  }
	    \item{seq}{path to mapped sample reads, formatted as either 'bed', 'tagAlign', or 'bowtie'}
	    \item{align}{path to directory containing alignability files for each chromsome, obtained from alignAdjust (check how these files are generated) or downloadd from our respository.  Alignability information is specific to the uniquness threshold one used to initially filter their mapped reads and the length of the sequence tags used}
	    \item{input}{path to mapped input reads, formatted as either 'bed', 'tagAlign', or 'bowtie'.  If left blank, then defaults to 'none'}
	    \item{extension}{average length of fragments in fragment library used, typically around 200}
	    \item{twobit}{path to current build of human genome in .2bit format}
	    \item{winSize}{window size, default=500bp}
	    \item{filetype}{format of mapped sample and input reads.  'bed' is files in the standard .bed format, 'tagAlign' signifies those in .taf format, and 'bowtie' signifies mapped reads directly outputted from bowtie.  Default is 'bowtie'}
	    \item{offset}{OPTIONAL. bp offset, default is none.  If one's winSize is 500 and would like 4 equally spaced offsets, then specifying offset=125 would achieve this}
	    \item{cnvWinSize}{OPTIONAL. Size of windows used to calculate CNV activity in sample, default size is 100000}
   	    \item{cnvOffset}{OPTIONAL.  Offset for CNV windows, typically 2500bp is suffcient, although default is no offsets}
   	    \item{outdir}{OPTIONAL.  Path to output directory if desired, otherwise places built files in same directory as the sample reads (outdir="default".}
	  }
          \seealso{
            \code{\link{save}}.
          }
          \examples{
#Chip-Seq with recommended parameters
buildwindowdata(
   filelist='data/ctcf.list', 
   seq='data/ctcfGM12878rep3chr22.taf' , 
   align='align1/', 
   input='data/inputGM12878rep3chr22.taf', 
   twoBit='hg18.2bit', 
   winSize=500, 
   offset=125, 
   cnvWinSize=100000, 
   cnvOffset=2500, 
   filetype='tagAlign',    
   extension=200)

#FAIRE-seq with recommended parameters, no input.  
buildwindowdata(
   filelist='data/faire.list', 
   seq='data/faireGM12878rep1chr22.taf', 
   align='align4/', 
   input='none', 
   twoBit='hg18.2bit', 
   winSize=250, 
   offset=50, 
   cnvWinSize=100000, 
   cnvOffset=2500, 
   filetype='tagAlign',    
   extension=134)


   
          }
          \keyword{file}
