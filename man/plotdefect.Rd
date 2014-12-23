\name{plotdefect}
\alias{plotdefect}

\title{
Mass defect vs. m/z scatterplot of HRMS peaks, with specific m/z isotope differences highlighted.
}

\description{
Mass defect vs. m/z scatterplot of HRMS peaks, with specific m/z isotope differences highlighted.
}

\usage{
plotdefect(pattern, elements = c("Br"))
}

\arguments{
  \item{pattern}{
	List of type pattern produced by \code{\link[nontarget]{pattern.search}}.
	}
  \item{elements}{
	Character string of an element for which the isotope m/z differences between two peaks detected by \code{\link[nontarget]{pattern.search}} should be highlighted (red).
}
}
\author{
Martin Loos
}
\note{
Here, mass defect is defined as the difference of m/z to the nearest integer from rounding.
\code{rm.comp = FALSE} and \code{rm.noncomp = FALSE} leads to no selection and thus no exclusion of anything.
}
\seealso{
	\code{\link[nontarget]{pattern.search}}
}
\examples{
\donttest{
data(peaklist);

peaklist<-rm.sat(peaklist,dmz=0.3,drt=0.1,intrat=0.015,spar=0.8,corcut=-1000,plotit=TRUE);

peaklist<-peaklist[peaklist[,4],1:3];

data(isotopes);

iso<-make.isos(isotopes,
	use_isotopes=c("13C","15N","34S","37Cl","81Br","41K","13C","15N","34S","37Cl","81Br","41K"),
	use_charges=c(1,1,1,1,1,1,2,2,2,2,2,2))
	
pattern<-pattern.search(
  peaklist,
  iso,
  cutint=10000,
  rttol=c(-0.05,0.05),
  mztol=2,
  mzfrac=0.1,
  ppm=TRUE,
  inttol=0.2,
  rules=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
  deter=FALSE,
  entry=50
);

plotdefect(pattern,elements=c("N"));
plotdefect(pattern,elements=c("Cl"));
plotdefect(pattern,elements=c("Br"));
plotdefect(pattern,elements=c("S"));
plotdefect(pattern,elements=c("C"));
plotdefect(pattern,elements=c("K"));
# P has only one isotope, hence:
# plotdefect(pattern,elements=c("P"));
}
}

