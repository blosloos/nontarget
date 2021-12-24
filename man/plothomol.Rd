\name{plothomol}
\alias{plothomol}
\title{
Marks homologue series peaks in a scatterplot of retention time (RT) vs. m/z
}
\description{
Given results from \code{\link[nontarget]{homol.search}}, a scatterplot of peaks within m/z and RT is generated with homologue series marked.
Herein, homologue series receive a color code based on the mean m/z differences between adjacent peaks of a series; these differences 
are rounded up to the second digit.
}
\usage{
plothomol(homol, xlim = FALSE, ylim = FALSE,plotlegend=TRUE,plotdefect=FALSE)
}

\arguments{
  \item{homol}{
	List of type homol produed by \code{\link[nontarget]{homol.search}}.
	}
  \item{xlim}{
	\code{xlim=c(upper bound,lower bound)}, default = FALSE.	
	}
  \item{ylim}{
	\code{ylim=c(upper bound,lower bound)}, default = FALSE.
	}
  \item{plotlegend}{
	Should a listing of m/z differences within homologue series and the concommittant color codes been added to the plot? If not, set to FALSE.
	}
   \item{plotdefect}{
	Plot the mass defect instead of the m/z value.
	}
}
\author{
Martin Loos
}

\seealso{
	\code{\link[nontarget]{homol.search}}
}
\examples{
\donttest{
data(peaklist);
data(isotopes)
homol<-homol.search(
	peaklist,
	isotopes,	
	elements=c("C","H","O"),
	use_C=TRUE,
	minmz=5,
	maxmz=120,
	minrt=2,
	maxrt=2,
	ppm=TRUE,
	mztol=3.5,
    rttol=0.5,
	minlength=3,
	mzfilter=FALSE,
	vec_size=3E6,
	spar=.45,
	R2=.98,
	plotit=FALSE
)
plothomol(homol,xlim=FALSE,ylim=FALSE,plotlegend=FALSE,plotdefect=FALSE);
}
}

