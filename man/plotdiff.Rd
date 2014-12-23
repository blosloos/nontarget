\name{plotdiff}
\alias{plotdiff}
\title{
Filtering important m/z differences among peaks of a HRMS data set.
}
\description{
Produce a vector and histogram of m/z differences among peaks in a HRMS data set. Frequent m/z differences may be relatable to isotope patterns and the presence of different adducts.
}
\usage{
plotdiff(peaklist, histbreaks = 10000, rttol = c(0, 0), mztol = c(0, 100), plotit = TRUE)
}
\arguments{
  \item{peaklist}{
	Dataframe of HRMS peaks with three numeric columns for (a) m/z, (b) intensity and (c) retention time, such as \code{\link[nontarget]{peaklist}}.
	}
  \item{histbreaks}{
	Number of histogram breaks.
	}
  \item{rttol}{
	Window (upper and lower difference bound, relative to the one peak) of retentiom time (RT) differences of peaks to the one peak screened from, see details and note.
	Units as given in column 3 of peaklist argument, e.g. [min].
	}
  \item{mztol}{
	Window (upper and lower difference bound, relative to the one peak) of m/z differences [u] of peaks to the one peak screened from, see details. 
	}
  \item{plotit}{
	Should histogram be plotted? If \code{FALSE}, \code{\link[nontarget]{plotdiff}} will only return a vector of type \code{diff}, see value. 
	}
}
\details{
	For each one peak in the dataset, \code{\link[nontarget]{plotdiff}} screens for other peaks within arguments \code{rttol} and \code{mztol}, saves their m/z difference to the m/z value 
	of the one peak and, over all one peaks, finally generates a histogram of all these m/z differences. Thus, and depending on the resolution set by argument \code{histbreaks}, 
	frequent m/z differences can be visualized.
}
\note{
	Argument \code{rttol} can e.g. be used to only include m/z differences of peaks with a higher RT relative to that of the one peak (as is the case in homologue series). 
	For example, let one peak have RT=12 min. Using \code{rttol=c(1,4)}, only m/z differences with peaks having a RT in between 13 and 14 min will then be screened for this one peak. 
	Akin for argument \code{mztol}.
}
\value{
	Vector \code{diffs} of m/z differences. Can serve as input to \code{\link[nontarget]{deter.iso}}.
}
\author{
Martin Loos
}
\seealso{
	\code{\link[nontarget]{peaklist}}
}
\examples{
\donttest{
data(peaklist)
diffs<-plotdiff(peaklist, histbreaks = 10000, rttol = c(0, 0), mztol = c(0, 100), plotit = TRUE)
}
}
