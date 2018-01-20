\name{summary.zinba}
\alias{summary.zinba}
\title{Printing summary of the ZINBA peak finding software}
\description{
    'run.zinba' runs the peak finder using the following parameters
}
\usage{
 run.zinba()

 buildwindowdata()

 getsigwindows()

 basecountimport()

 peakbound()
}


\arguments{
  \item{x}{an object from the 'zicounts' or 'zicensor' fit}
  \item{object}{an object from the 'zicounts' or 'zicensor' fit}
  \item{digits}{number of significant decimal points.}  
  \item{...}{other options 'summary' arguments.}  

}

\references{
Bohning, D., E. Dietz, P. Schlattman, L. Mendonca, and U. Kirchner (1999). 
	The zero-inflated Poisson model and the decayed, missing and filled 
	teeth index in dental epidemiology. Journal of the Royal Statistical Society,
	Series A 162, 195--209.
	
Carrivick, P. J. W., A. H. Lee, and K. K. W. Yau (2003). Zero-inflated Poisson 
	to evaluate occupational safety interventions.
	Safety Science 41, 53--63.

Lambert, D. (1992). Zero-inflated Poisson regression, with an application to defects in
	manufacturing. Technometrics 34, 1--14.

Lewsey, J. D. and W. M. Thomson (2004). The utility of the zero-inflated Poisson and
	zero-inflated negative binomial models: a case study of cross-sectional and
	longitudinal DMF data examining the effect of socio-economic status. 
	Community Dentistry and Oral Epidemiology 32, 183--189.
	
Mwalili, S. M., Lesaffre, E. and Declerck, D. (2004). The interval censored zero-inflated
	negative binomial regression model: an application in caries research. submitted.

}

\seealso{
\code{\link{zicounts}}.
}

\examples{
## See help for ?zinba
}

\keyword{models}
\author{Naim Rashid <nrashid@bios.unc.edu>}
