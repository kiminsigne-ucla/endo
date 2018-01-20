run=function(util){
	###################################################################
	###################################################################
	if(util=="Generate Alignability"){
	genalign <- aDialog(items=list(
		mapdir=fileItem(attr=c(type="selectdir"), tooltip="Path to unpacked mappability directory", label="Mappability\nDirectory"),	
		outdir=fileItem(attr=c(type="selectdir"), tooltip="Path to alignability directory where output files will be placed",label="Output\nDirectory"), 
		athresh=numericItem(1, tooltip="Alignability threshold used to filter mapped reads, 1 corresponding that reads only mapping to one place in genome were used", label="Alignability\nThreshold"), 
		extension=numericItem(200, tooltip="Average length in fragment library used in preparing sample/input reads for sequencing", label="Extension"),  
		twoBitFile=fileItem("", attr=list(
        	                             filter=list(
					       "2bit" = list(patterns=c("*.2bit")),
        	                               "All files" = list(patterns=c("*"))
        	                               )), tooltip="Select .2bit file for analysis", label=".2bit file")
		), title="Generate Alignability Directory")

	                           
	view_genalign <- aContainer("mapdir", "outdir", "athresh", "extension", "twoBitFile")
	genalign$make_gui(gui_layout=view_genalign, visible=FALSE) 
	genalign$OK_handler <- function(.) {
		genalignvar=.$to_R()
		do.call("generateAlignability",list(mapdir=genalignvar$mapdir, outdir=genalignvar$outdir, athresh=genalignvar$athresh, extension=genalignvar$extension, twoBitFile=genalignvar$twoBitFile)) 		
		.$close_gui()
	}
		genalign$visible(TRUE)
	}

	if(util=="Generate SBC file (basecount)"){	
	basecount <- aDialog(items=list(
		seq=fileItem("", attr=list(
                                     filter=list("Bed"=list(
                                                   patterns=c("*.bed")
                                                   ),
				       "tagAlign" = list(patterns=c("*.taf","*.tagAlign")),
                                       "All files" = list(patterns=c("*"))
                                       )),tooltip="Select mapped sample reads file for analysis",  label="Sample Reads"),	
		basecountdir=fileItem(attr=c(type="selectdir"), tooltip="Directory where basecountfile will be placed", label="Output\nDirectory"),			
		filename=stringItem("name of SBC (basecount) file to be created", tooltip="Name of SBC (basecount) file, for example ctcf.basecount", label="SBC file"),
		twoBitFile=fileItem("", attr=list(
                                     filter=list(
				       "2bit" = list(patterns=c("*.2bit")),
                                       "All files" = list(patterns=c("*"))
                                       )), tooltip="Select .2bit file for analysis", label=".2bit file"),
		extension=numericItem(200, tooltip="Average length in fragment library used in preparing sample/input reads for sequencing", label="Extension"),  
		filetype=choiceItem("bed", values=c("bed","tagAlign", "bowtie"), editor_type="gcombobox",tooltip="Select format of mapped sample read file", label="File Type")
		), title="Generate SBC (basecount) File")
	                           
	view_basecount <- aContainer("seq","filetype","basecountdir","filename","twoBitFile","extension")
	basecount$make_gui(gui_layout=view_basecount, visible=FALSE) 
	basecount$OK_handler <- function(.) {
		basecountvar=.$to_R()
		do.call("basealigncount",list(inputfile=basecountvar$seq,outputfile=paste(basecountvar$basecount, "/",basecountvar$filename,sep="") ,twoBitFile=basecountvar$twoBitFile,extension=basecountvar$extension,filetype=basecountvar$filetype))
		.$close_gui()
	}
		basecount$visible(TRUE)
	}
	
	if(util=="twobittofa"){	
	twobittofa <- aDialog(items=list(
		twoBitFile=fileItem("", attr=list(
                                     filter=list(
				       "2bit" = list(patterns=c("*.2bit")),
                                       "All files" = list(patterns=c("*"))
                                       )), tooltip="Select .2bit file for Conversion to fasta format",label=".2bit file"),
		chrm=stringItem("all", tooltip="Name of chromosome where sequence will be converted from, or all if convering an entire directory of fasta files", label="chromosome"),
		outdir=fileItem(attr=c(type="selectdir"), tooltip="Path to directory where converted fasta files will be placed",label="Output\nDirectory"), 
		start=numericItem(0, tooltip="Base Start Position (if chromosome != all) of sequence to be extracted to fasta format",label="Sequence Start"),  
		end=numericItem(0, tooltip="Base Stop Position (if chromosome != all) of sequence to be extracted to fasta format",label="Sequence End"),
		gcSeq=stringItem("path to output file", tooltip="Path to file where extracted fasta sequence will be placed", label="chromosome")
		), title="Convert .2bit to fasta format")
	                           
	view_twobittofa <- aContainer("twoBitFile","chrm","outdir",
					aContext("start", context=twobittofa, 
						enabled_when=function(.) {
						val <- .$to_R()$chrm
						val !="all"}),
					aContext("end", context=twobittofa, 
								enabled_when=function(.) {
  								val <- .$to_R()$chrm
  								val !="all"
							}),
					aContext("gcSeq", context=twobittofa, 
								enabled_when=function(.) {
  								val <- .$to_R()$chrm
  								val !="all"
							})
					)
	twobittofa$make_gui(gui_layout=view_twobittofa, visible=FALSE) 
	twobittofa$OK_handler <- function(.) {
		twobittofavar=.$to_R()
		do.call("twobittofa",list(twoBitFile=twobittofavar$twoBitFile,chrm=twobittofavar$chrm,outdir=twobittofavar$outdir, start=twobittofavar$start,end=twobittofavar$end,gcSeq=twobittofavar$gcSeq))
		.$close_gui()
	}
		twobittofa$visible(TRUE)
	}

	if(util=="fatotwobit"){	
		fatotwobit <- aDialog(items=list(
			outFile=stringItem("/home/data/genome.2bit", tooltip="if a directory is not selected, path to file where resulting .2bit file will be outputted", label="Path to\nOutput File"),
			fadir=fileItem(" ",attr=c(type="selectdir"), tooltip="Path to directory where fasta files to be converted are stored", label="fasta directory"), 
			faFile=stringItem(" ", tooltip="path to specific fasta file if not specifying a directory", label="chromosome")
			), title="Convert fasta file to to .2bit format")
		                        
		view_fatotwobit <- aContainer("fadir",
						aContext("faFile", context=fatotwobit, 
							enabled_when=function(.) {
							val <- .$to_R()$fadir
							val ==" "}),
						"outFile")
		fatotwobit$make_gui(gui_layout=view_fatotwobit, visible=FALSE) 
		fatotwobit$OK_handler <- function(.) {
			fatotwobitvar=.$to_R()
			do.call("fatotwobit",list(fadir=fatotwobitvar$fadir,faFile=fatotwobitvar$faFile,outFile=fatotwobitvar$outFile))
			.$close_gui()
		}
			fatotwobit$visible(TRUE)
		}
	if(util=="getwindowcounts"){	
		getwindowcounts <- aDialog(items=list(
		 	seq=fileItem("", attr=list(
                 	                    filter=list("Bed"=list(
                 	                                  patterns=c("*.bed")
                 	                                  ),
					       "tagAlign" = list(patterns=c("*.taf","*.tagAlign")),
                 	                      "All files" = list(patterns=c("*"))
                 	                      )),tooltip="Select mapped sample reads file for analysis", label="Sample Reads"),
		 	filetype=choiceItem("bed", values=c("bed","tagAlign", "bowtie"), editor_type="gcombobox",tooltip="Select format of sample and input mapped read file",label="File Type"),
	         	twoBit=fileItem("", attr=list(
                 	                    filter=list(
					       "2bit" = list(patterns=c("*.2bit")),
                 	                      "All files" = list(patterns=c("*"))
                 	                      )), tooltip="Select .2bit file corresponding on build of genome reads were mapped to",label=".2bit file"),	
                 winSize=integerItem(500, tooltip="Window Size for data, smaller is better for low signal to noise, high background data", label="Window Size"),
                 offset=integerItem(0, tooltip="BP offset to shift windows, window size must be a multiple of this number.  Yields better sensitivity", label="Window Offset"),
                 extension=integerItem(200, tooltip="average fragment length in sample library used for sequencing reads, usually 150bp-200bp", label="Extension")
		),title="Get Window Read Counts for data")
                 
		                           
		view_getwindowcounts <- aContainer("seq","filetype","twoBit","winSize","offset", "extension")
		getwindowcounts$make_gui(gui_layout=view_getwindowcounts, visible=FALSE) 
		getwindowcounts$OK_handler <- function(.) {
			getwindowcountsvar=.$to_R()
			do.call("getwindowcounts",list(seq=getwindowcountsvar$seq,filetype=getwindowcountsvar$filetype,twoBit=getwindowcountsvar$twoBit,winSize=getwindowcountsvar$winSize, offset=getwindowcountsvar$offset, extension=getwindowcountsvar$extension))
			.$close_gui()
		}
			getwindowcounts$visible(TRUE)
		}
if(util=="Signal to Noise Analysis"){	
		signalnoise <- aDialog(items=list(
		 	inputfile=fileItem("", attr=list(
                 	                    filter=list("Bed"=list(
                 	                                  patterns=c("*.bed")
                 	                                  ),
					       "basecount" = list(patterns=c("*.basecount")),
                 	                      "All files" = list(patterns=c("*"))
                 	                      )),tooltip="Select basecount file for analysis", label="SBC (basecount) File"),
		 	twoBit=fileItem("", attr=list(
                 	                    filter=list(
					       "2bit" = list(patterns=c("*.2bit")),
                 	                      "All files" = list(patterns=c("*"))
                 	                      )), tooltip="Select .2bit file corresponding on build of genome reads were mapped to",label=".2bit file"),	
                winSize=integerItem(100000, tooltip="Window Size for analysis, larger, 100kb suggested", label="Window Size")), title="Signal to Noise Analysis")
                
		view_signalnoise <- aContainer("inputfile","twoBit","winSize")
		signalnoise$make_gui(gui_layout=view_signalnoise, visible=FALSE) 
		signalnoise$OK_handler <- function(.) {
			signalnoisevar=.$to_R()
			do.call("signalnoise",list(inputfile=signalnoisevar$inputfile,twoBit=signalnoisevar$twoBit,winSize=signalnoisevar$winSize))
			.$close_gui()
		}
			signalnoise$visible(TRUE)
		}
	###################################################################
	###################################################################
	
	if(util=="Run Zinba"){
		dlg <- aDialog(items=list(
                 RunChoice=choiceItem("", values=c("Build Windows?","Run Window Analysis?", "Refine Merged Significant Regions?"), tooltip="The Zinba steps you want to run, either together or alone", editor_type="gcheckboxgroup", , label="    Run Options"),
		 printFullOut=choiceItem(as.integer(1), values=c(as.integer(1), as.integer(0)),tooltip="If 0, only prints out window coordinates, filter status, and window enrichment probability columns for .wins files from window analysis, rather than whole dataset",editor_type="gcombobox", label="   Print Option"),
		 method=choiceItem("mixture", values=c("mixture", "pscl"), tooltip="analysis method to use, either mixture regression model or pscl method (FAIRE)",editor_type="gcombobox", label="Analysis Method"),
		 outfile=stringItem("prefix for output files", tooltip="path and prefix to where .wins and .peaks are placed if analysis and refinement selected together, for example '/home/data/analysis_1'"),	
 		 numProc=integerItem(1, tooltip="How many processing cores to allocate, more processors means faster execution but more memory needed", label="Num Processors"),	 
		 twoBit=fileItem("", attr=list(
                                     filter=list(
				       "2bit" = list(patterns=c("*.2bit")),
                                       "All files" = list(patterns=c("*"))
                                       )), tooltip="Select .2bit file corresponding on build of genome reads were mapped to",label=".2bit file"),
                 seq=fileItem("", attr=list(
                                     filter=list("Bed"=list(
                                                   patterns=c("*.bed")
                                                   ),
				       "tagAlign" = list(patterns=c("*.taf","*.tagAlign")),
                                       "All files" = list(patterns=c("*"))
                                       )),tooltip="Select mapped sample reads file for analysis", label="Sample Reads"),
                 input=fileItem("none", attr=list(
                                     filter=list("Bed"=list(
                                                   patterns=c("*.bed")
                                                   ),
				       "tagAlign" = list(patterns=c("*.taf","*.tagAlign")),
                                       "All files" = list(patterns=c("*"))
                                       )),tooltip="Select mapped input control reads file (if available) for analysis",label="Input Reads"),
		 filetype=choiceItem("bed", values=c("bed","tagAlign", "bowtie"), editor_type="gcombobox",tooltip="Select format of sample and input mapped read file",label="File Type"),
		 align=fileItem(attr=c(type="selectdir", text="Path To Align Directory"), label="Align Directory"),	
                 createfilelist=stringItem("Path To filelist", tooltip="Path to filelist that will be created, holds locations of build datasets"),
		 winSize=integerItem(500, tooltip="Window Size for data, smaller is better for low signal to noise, high background data", label="Window Size"),
                 offset=integerItem(125, tooltip="BP offset to shift windows, window size must be a multiple of this number.  Yields better sensitivity", label="Window Offset"),
                 extension=integerItem(200, tooltip="average fragment length in sample library used for sequencing reads, usually 150bp-200bp", label="Extension"),
                 filelist=fileItem("Path existing file list", attr=list(
                                     filter=list(
				       "file list" = list(patterns=c("*.list")),
                                       "All files" = list(patterns=c("*"))
                                       )), label="Select Filelist"),
		 formula=stringItem("", tooltip="formula for modeling background: some combination of gcPerc, align_perc, exp_cnvwin_log, and input_count", label="     BG Formula"),
		 formulaE=stringItem("exp_count~1", tooltip="formula for modeling enriched region: some combination of, gcPerc, align_perc, exp_cnvwin_log, and input_count", label="   Enrich Formula"),
		 threshold=numericItem(.01, tooltip="FDR threshold for the pscl method only, lower is more stringent", label="PSCL Thresh"),
	         peakconfidence=numericItem(.8, tooltip="Posterior window encrichment probability threshold for mixture method only, higher is more stringent", label="Mixture Thresh"),
	         basecountfile=fileItem("", attr=list(
                                     filter=list(
				       "basecount" = list(patterns=c("*.basecount")),
                                       "All files" = list(patterns=c("*"))
                                       )), tooltip="Path to SBC (basecountfile) generated from Generate SBC (Basecount) File()"),
		 winlist=fileItem("Path to.winlist from previous window analysis", attr=list(
                                     filter=list(
				       "winlist" = list(patterns=c("*.winlist")),
                                       "All files" = list(patterns=c("*"))
                                       )), label="        winlist"),
		 peaksfilename=stringItem("Path To file for refined peaks", tooltip="Path To file for refined peaks"),
		 pWinSize=integerItem(200, tooltip='sliding window size to detect local maximums, smaller values are more sensitive but may yield insignificant peaks'),
		 pquant=numericItem(1, tooltip='quantile threshold for basecounts, 1 means taking only the max of the entire merged region'),
		 cnvWinSize=integerItem(100000, tooltip="Window size for estimating cnv activity"),
                 cnvOffset=integerItem(2500, tooltip="BP offset to shift CNV windows, CNV window size must be a multiple of this number"),
		 tol=numericItem(.00001, tooltip='mixture regression convergence relative tolerance level'),
		 initmethod=choiceItem("count", values=c("count","quantile", "pscl"),editor_type="gcombobox") ,
 		cleanup=trueFalseItem(FALSE, tooltip=paste("If true, then intermediate files are deleted once zinba is complete"))         
                 

  ),
title="Zinba Custom GUI",
help_string="Select a file, then adjust parameters.")

view <- aNotebook(
                  aNotebookPage( aFrame(
					aFrame(aContainer("RunChoice" ,"printFullOut","method",aContext("threshold", context=dlg, enabled_when=function(.) {
  								val <- .$to_R()$method
  								val =="pscl"
							}),aContext("peakconfidence", context=dlg, 
								enabled_when=function(.) {
  								val <- .$to_R()$method
  								val =="mixture"
							}),"outfile","numProc","twoBit", "cleanup"),label='General Options'
					), 

                                        aFrame(aContainer("seq","input","filetype","align","createfilelist","winSize","offset", "extension"),label='Build Window Data',
							enabled_when=function(.) {
  								val <- .$to_R()$RunChoice
  								sum(val=="Build Windows?")==1
							} 
					)
					,label='', horizontal=TRUE
				),

				aFrame(
					aFrame(aContainer(
							aContext("filelist", context=dlg, 
								enabled_when=function(.) {
  								val <- .$to_R()$RunChoice
  								sum(val =="Build Windows?")==0
							})									
							,"formula", 
							aContext("formulaE", context=dlg, 
								enabled_when=function(.) {
  								val <- .$to_R()$method
  								val =="mixture"
							}), "cnvWinSize", "cnvOffset","tol", "initmethod"),
							label='Window Analysis Options',
							enabled_when=function(.) {
  								val <- .$to_R()$RunChoice
  								sum(val=="Run Window Analysis?")==1
							}
					), 
					aFrame(aContainer("basecountfile",
							aContext("winlist", context=dlg, 
								enabled_when=function(.) {
  								val <- .$to_R()$RunChoice
  								sum(val =="Build Windows?")==0 & sum(val=="Run Window Analysis?")==0
							}),
							aContext("peaksfilename", context=dlg, 
								enabled_when=function(.) {
  								val <- .$to_R()$RunChoice
  								sum(val =="Build Windows?")==0 & sum(val=="Run Window Analysis?")==0
							})
							,"pWinSize","pquant"),label='Peak Refinement Options', 
						enabled_when=function(.) {
  								val <- .$to_R()$RunChoice
  								sum(val=="Refine Merged Significant Regions?")==1
							}
					), 
					label='', horizontal=TRUE
				),
			label="Run Zinba"
			)       
		)
		
	dlg$OK_handler <- function(.) {
		var=.$to_R()
		#set up exectution options
		if(sum(var$RunChoice=="Build Windows?")==1){
			listvar=var$createfilelist
			build=1
		}else{
			listvar=var$filelist
			build=0
		}
		if(sum(var$RunChoice=="Run Window Analysis?")==1){
			analysis=1
		}else{
			analysis=0
		}
		if(sum(var$RunChoice=="Refine Merged Significant Regions?")==1){
			refine=1
		}else{
			refine=0
		}
		
		#run conditional execution
		if(build==1 & analysis==0 & refine==0){
			do.call("buildwindowdata",list(seq=var$seq,align=var$align,input=var$input,twoBit=var$twoBit,winSize=var$winSize,offset=var$offset,cnvWinSize=var$cnvWinSize,cnvOffset=var	$cnvOffset,filelist=listvar,filetype=var$filetype,extension=var$extension))
		}else if(build==0 & analysis==0 & refine==1){
			do.call("getrefinedpeaks",list(winlist=var$winlist,basecountfile=var$basecountfile,bpout=paste(var$outfile,'.temp', sep=''),peakout=var$peaksfilename,twoBit=var$twoBit,winSize=500,pWinSize=var$pWinSize,pquant=var$pquant,printFullOut=var$printFullOut,peakconfidence=var$peakconfidence,threshold=var$threshold,method=var$method))
		}else{
			do.call("run.zinba",list(filelist=listvar,formula=as.formula(var$formula),formulaE=as.formula(var$formulaE), outfile=var$outfile, seq=var$seq,align=var$align,input=var$input,twoBit=var$twoBit,winSize=var$winSize,offset=var$offset,cnvWinSize=var$cnvWinSize,cnvOffset=var$cnvOffset, basecountfile=var$basecountfile, threshold=var$threshold,peakconfidence=var$peakconfidence,tol=var$tol,numProc=var$numProc, buildwin=build, pWinSize=var$pWinSize,pquant=var$pquant,refinepeaks=refine, printFullOut=var$printFullOut,method=var$method,initmethod=var$initmethod, diff=0, filetype=var$filetype,extension=var$extension, cleanup=var$cleanup ))
		}
 	}
	dlg$make_gui(gui_layout=view)
	}
	###################################################################
	###################################################################
if(util=="Run Zinba Preset"){
		dlg <- aDialog(items=list(
                 preset=choiceItem("High Signal-to-Noise", values=c("High Signal-to-Noise","Moderate STN/Broad-Sharp Mix","Low STN/Broad"), tooltip="Examples:High STN = CTCF ChIP-seq; Moderate STN/Broad-Sharp Mix = PolII ChIP-seq; Low STN = FAIRE-seq/Histone-seq", editor_type="gcombobox", , label="    Data Type"),
		 method=choiceItem("mixture", values=c("mixture", "pscl"), tooltip="analysis method to use, either mixture regression model or pscl method (FAIRE)",editor_type="gcombobox", label="Analysis Method"),
		 outfile=stringItem("prefix for output files", tooltip="path and prefix to where .wins and .peaks are placed if analysis and refinement selected together, for example '/home/data/analysis_1'"),	
 		 numProc=integerItem(1, tooltip="How many processing cores to allocate, more processors means faster execution but more memory needed", label="Num Processors"),	 
		 twoBit=fileItem("", attr=list(
                                     filter=list(
				       "2bit" = list(patterns=c("*.2bit")),
                                       "All files" = list(patterns=c("*"))
                                       )), tooltip="Select .2bit file corresponding on build of genome reads were mapped to",label=".2bit file"),
                 seq=fileItem("", attr=list(
                                     filter=list("Bed"=list(
                                                   patterns=c("*.bed")
                                                   ),
				       "tagAlign" = list(patterns=c("*.taf","*.tagAlign")),
                                       "All files" = list(patterns=c("*"))
                                       )),tooltip="Select mapped sample reads file for analysis", label="Sample Reads"),
                 input=fileItem("none", attr=list(
                                     filter=list("Bed"=list(
                                                   patterns=c("*.bed")
                                                   ),
				       "tagAlign" = list(patterns=c("*.taf","*.tagAlign")),
                                       "All files" = list(patterns=c("*"))
                                       )),tooltip="Select mapped input control reads file (if available) for analysis",label="Input Reads"),
		 filetype=choiceItem("bed", values=c("bed","tagAlign", "bowtie"), editor_type="gcombobox",tooltip="Select format of sample and input mapped read file",label="File Type"),
		 align=fileItem(attr=c(type="selectdir", text="Path To Align Directory"), label="Align Directory"),	
                 extension=integerItem(200, tooltip="average fragment length in sample library used for sequencing reads, usually 150bp-200bp", label="Extension"),
                 basecountfile=fileItem("", attr=list(
                                     filter=list(
				       "basecount" = list(patterns=c("*.basecount")),
                                       "All files" = list(patterns=c("*"))
                                       )), tooltip="Path to SBC (basecountfile) generated from Generate SBC (Basecount) File()"),
		cleanup=trueFalseItem(FALSE, tooltip=paste("If true, then intermediate files are deleted once zinba is complete"))
  ),
title="Run Zinba Pipeline With Data Presets",
help_string="Select a file, then adjust parameters.")

view <- aFrame(aContainer("preset" ,"method","outfile","numProc","twoBit","seq",aContext("input",context=dlg, 
											enabled_when=function(.) {
				  								val <- .$to_R()$preset
  												sum(val ==c("High Signal-to-Noise","Moderate STN/Broad-Sharp Mix"))>0
											}),
			"filetype","align","extension","basecountfile", "cleanup"), 
			label='', horizontal=TRUE)
		       
		
		
dlg$OK_handler <- function(.) {
	var=.$to_R()		
	if(var$preset=="High Signal-to-Noise"){
		printFullOut=1;
		 winSize=500;
                 offset=125;
                 formula="exp_count~input_count";
		 formulaE="exp_count~1";
		 threshold=.01;
	         peakconfidence=.99;
	         pWinSize=200;
		 pquant=1;
		 cnvWinSize=100000;
                 cnvOffset=2500;
		 tol=.00001;
		 initmethod="count"
	}else if(var$preset=="Moderate STN/Broad-Sharp Mix"){
		printFullOut=1;
		 winSize=1000;
                 offset=250;
                 formula="exp_count~input_count";
		 formulaE="exp_count~1";
		 threshold=.01;
	         peakconfidence=.99;
	         pWinSize=200;
		 pquant=1;
		 cnvWinSize=100000;
                 cnvOffset=2500;
		 tol=.00001;
		 initmethod="count"
	}else if(var$preset=="Low STN/Broad"){
		printFullOut=1;
		 winSize=250;
                 offset=125;
                 formula="exp_count~gcPerc+align_perc+exp_cnvwin_log";
		 formulaE="exp_count~1";
		 threshold=.01;
	         peakconfidence=.8;
	         pWinSize=200;
		 pquant=1;
		 cnvWinSize=100000;
                 cnvOffset=2500;
		 tol=.00001;
		 initmethod="count"		
		}
		finalrunlist=list(filelist=NULL,formula=as.formula(formula),formulaE=as.formula(formulaE), outfile=var$outfile, seq=var$seq,align=var$align,input=var$input,twoBit=var$twoBit,winSize=winSize,offset=offset,cnvWinSize=cnvWinSize,cnvOffset=cnvOffset, basecountfile=var$basecountfile, threshold=threshold,peakconfidence=peakconfidence,tol=tol,numProc=var$numProc, buildwin=1, pWinSize=pWinSize,pquant=pquant,refinepeaks=1, printFullOut=printFullOut,method=var$method,initmethod=initmethod, diff=0, filetype=var$filetype,extension=var$extension, cleanup=var$cleanup)
		print(finalrunlist)
		do.call("run.zinba",finalrunlist)
		
 	}
	dlg$make_gui(gui_layout=view)
	}
}


zinbagui=function(){
require(traitr)
require(gWidgets)
require(zinba)
options(guiToolkit="RGtk2")
eval(parse(text="statusold=c('Select a utility','Yes','Yes', '?')"), envir=.GlobalEnv)
dlg <- aDialog(items=list(
	util=choiceItem("Select a utility", values=c("Select a utility","Generate Alignability Directory", "Generate SBC (Basecount) File", "Convert .2bit to .fa","Convert .fa to .2bit", "Get Window Counts","Signal to Noise Analysis"), tooltip="Various utilities to compute files needed for ZINBA",editor_type="gcombobox", label="Utilities"),
	align=choiceItem("Yes", values=c("Yes", "No"), tooltip="If no, then will open a dialog to generate your alignability directory for you",editor_type="gcombobox", label="Have you generated your alignability directory?"),
	basecount=choiceItem("Yes", values=c("Yes", "No"), tooltip="If no, then will open a dialog to generate your SBC (basecount) file",editor_type="gcombobox",label="Have you generated your SBC (basecount) track?"),
	runzinba=choiceItem("?", values=c("?","Run Zinba Data Preset","Run Zinba Custom Parameters"), tooltip="Open Zinba Dialog Window",editor_type="gcombobox", label="Zinba Options")
	),
	title=" Zinba GUI ",
	model_value_changed=function(.) {
                 l <- .$to_R()
		 if(l$util=="Convert .2bit to .fa" & l$util!=statusold[1]){
			do.call("run", list("twobittofa"))
			eval(parse(text=paste("statusold[1]='",l$util,"'", sep="")), envir=.GlobalEnv)
		}else if(l$util=="Generate Alignability Directory" & l$util!=statusold[1]){
			do.call("run", list("Generate Alignability"))
			eval(parse(text=paste("statusold[1]='",l$util,"'", sep="")), envir=.GlobalEnv)
		}else if(l$util=="Generate SBC (Basecount) File" & l$util!=statusold[1]){
			do.call("run", list("Generate SBC file (basecount)"))
			eval(parse(text=paste("statusold[1]='",l$util,"'", sep="")), envir=.GlobalEnv)
		}else if(l$util=="Convert .fa to .2bit" & l$util!=statusold[1]){
			do.call("run", list("fatotwobit"))
			eval(parse(text=paste("statusold[1]='",l$util,"'", sep="")), envir=.GlobalEnv)
		}else if(l$util=="Get Window Counts" & l$util!=statusold[1]){
			do.call("run", list("getwindowcounts"))
			eval(parse(text=paste("statusold[1]='",l$util,"'", sep="")), envir=.GlobalEnv)
		}else if(l$util=="Signal to Noise Analysis" & l$util!=statusold[1]){
			do.call("run", list("Signal to Noise Analysis"))
			eval(parse(text=paste("statusold[1]='",l$util,"'", sep="")), envir=.GlobalEnv)
		}else if(l$util=="Select a utility" & l$util!=statusold[1]){
			eval(parse(text=paste("statusold[1]='",l$util,"'", sep="")), envir=.GlobalEnv)
		}

		if(l$align=="No" & l$align!=statusold[2]){
			do.call("run", list("Generate Alignability"))
			eval(parse(text=paste("statusold[2]='",l$align,"'", sep="")), envir=.GlobalEnv)
		}			
		if(l$basecount=="No" & l$basecount!=statusold[3]){
			do.call("run", list("Generate SBC file (basecount)"))	
			eval(parse(text=paste("statusold[3]='",l$basecount,"'", sep="")), envir=.GlobalEnv)	
	 	}
		if(l$runzinba =="Run Zinba Custom Parameters" & l$runzinba!=statusold[4]){
			do.call("run", list("Run Zinba"))
			eval(parse(text=paste("statusold[4]='",l$runzinba,"'", sep="")), envir=.GlobalEnv)
		 }else if(l$runzinba =="Run Zinba Data Preset"& l$runzinba!=statusold[4]) 
			do.call("run", list("Run Zinba Preset"))
			eval(parse(text=paste("statusold[4]='",l$runzinba,"'", sep="")), envir=.GlobalEnv)				
               }
)

view <- aContainer(aFrame("util", label="Zinba Utilities"),aFrame("align", label="File Check"), aFrame("basecount", label="File Check"), aFrame("runzinba", label="Choose How to Run Zinba"))
dlg$make_gui(gui_layout=view)

}
