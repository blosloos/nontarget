\name{make.isos}
\alias{make.isos}

\title{
Deriving list of m/z isotope differences for input into \code{\link[nontarget]{pattern.search}}.
}

\description{
Deriving list of m/z isotope differences for input into \code{\link[nontarget]{pattern.search}}.
}

\usage{
make.isos(isotopes,
use_isotopes=c("13C","15N","34S","37Cl","81Br","41K","13C",
"15N","34S","37Cl","81Br","41K"),
use_charges=c(1,1,1,1,1,1,2,2,2,2,2,2))
}

\arguments{
  \item{isotopes}{
	Dataframe with isotopes, from dependency \code{enviPat}
}
  \item{use_isotopes}{
	Character string of non-monoisotopic isotopes for isotopologue search. 
}
  \item{use_charges}{
	Numeric vector with length equal to that of \code{use_isotopes},. Specifies the charge z for the isotopologue search of each isotope.
}
}


\value{

List of type iso with 5 entries

\item{iso[[1]]}{\code{list of isotopes}.}
\item{iso[[2]]}{\code{list of isotope masses}.}
\item{iso[[3]]}{\code{charges}.}
\item{iso[[4]]}{\code{number of isotope m/z}.}
\item{iso[[5]]}{\code{elements}.}

}

\author{
Martin Loos
}

\seealso{
	\code{\link[nontarget]{pattern.search}}
}
\examples{
data(isotopes)
iso<-make.isos(isotopes,
	use_isotopes=c("13C","15N","34S","37Cl","81Br","41K","13C","15N","34S","37Cl","81Br","41K"),
	use_charges=c(1,1,1,1,1,1,2,2,2,2,2,2))
}

