\name{plotall}
\alias{plotall}

\title{
RT vs. m/z scatterplot marking isotope pattern and adduct group peaks.
}
\usage{
plotall(pattern, adduct)
}

\description{
RT vs. m/z scatterplot marking isotope pattern and adduct group peaks.
}

\arguments{
  \item{pattern}{
	List of type pattern produced by \code{\link[nontarget]{pattern.search}}.
	}
  \item{adduct}{
	List of type adduct produced by \code{\link[nontarget]{adduct.search}}.
	}
}
\author{
Martin Loos
}
\seealso{
	\code{\link[nontarget]{pattern.search}}
	\code{\link[nontarget]{adduct.search}}
}
\examples{
\donttest{

data(peaklist);
data(adducts);
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
adduct<-adduct.search(
  peaklist,
  adducts,
  rttol=0.05,
  mztol=3,
  ppm=TRUE,
  use_adducts=c("M+K","M+H","M+Na","M+NH4"),
  ion_mode="positive"
);
plotall(pattern, adduct)

}
}
