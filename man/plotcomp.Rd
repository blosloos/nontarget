\name{plotcomp}
\alias{plotcomp}

\title{
Plot and print isotope and adduct relations among peaks of a single component
}
\description{
Plot and print isotope and adduct relations among peaks of a single component. Also lists all other peaks of the 
data set within tolerance ranges of m/z and retention time (RT). 
}
\usage{
plotcomp(comp, compoID, peakID = FALSE)
}

\arguments{
  \item{comp}{
	List of type \code{comp} produced by function \code{\link[nontarget]{combine}}.
	}
  \item{compoID}{
	ID of component to be plotted. For description of component IDs see \code{\link[nontarget]{combine}}, note section. 
	Use with argument \code{peakID=FALSE}.
}
  \item{peakID}{
	ID of a peak in a component; selects the component containing the peak with this ID. For description of peak IDs see note section. 
	Use with argument \code{compoID=FALSE}.
	}
}
\details{
The upper plot panel provides a circular plot of peak relations, with m/z increasing clockwise starting from noon. Herein, peaks are 
represented by their peak IDs; numbers in brackets give decreasing peak intensity ranks over all peaks in the shown component.
Adduct relations are symbolized by red lines and isotope relations by blue arrows. Thin instead of thick lines stand for interfering peaks.
In addition, all relations, other peaks within range and homologue series information are printed as value of \code{\link[nontarget]{plotcomp}}

The lower panel barplot shows intensities vs. m/z of both the peaks in the component (bold) and the peaks within tolerance ranges of m/z and RT (grey), 
defined by arguments \code{mztol} and \code{rttol} of \code{\link[nontarget]{pattern.search}} and \code{\link[nontarget]{adduct.search}}.

}
\author{
Martin Loos
}
\note{
Input peaklist is internally sorted and saved in the lists returned by (a) increasing retention time and (b) m/z by all \code{\link[nontarget]{pattern.search}}, \code{\link[nontarget]{adduct.search}} 
and \code{\link[nontarget]{homol.search}}. Peak IDs refer to this very order - in contrast to group IDs. Different IDs exist for adduct groups, isotope pattern groups, grouped homologue series (HS) peaks 
and homologue series cluster. Moreover, and at the highest level, IDs exist for the individual components (see note section of \code{\link[nontarget]{combine}}).

}

\seealso{
	\code{\link[nontarget]{combine}}
}
