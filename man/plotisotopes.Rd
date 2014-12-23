\name{plotisotopes}
\alias{plotisotopes}
\title{
Plot of isotope counts over isotope pattern groups or components.
}
\description{
Plots and prints counts of m/z isotope differences detected either by \code{\link[nontarget]{pattern.search}} or by \code{\link[nontarget]{combine}}.
}
\usage{
plotisotopes(input)
}
\arguments{
  \item{input}{
	Either list of type \code{pattern} produced by \code{\link[nontarget]{pattern.search}} or a list of type \code{comp} produced by \code{\link[nontarget]{combine}}
	}
}
\details{
The function allows to track the number of m/z isotope differences (a) over individual pairs of peaks and (b) aggregated over isotope pattern groups (argument \code{pattern})
or (c) aggregated over components and (d) aggregated over components within small mass tolerance (argument \code{comp}).
The small mass tolerance is defined by the \code{massfrac} and \code{mztol} arguments of \code{\link[nontarget]{pattern.search}} and \code{\link[nontarget]{adduct.search}}.
}
\value{
Dataframe listing counts
}
\author{
Martin Loos
}
\seealso{
	\code{\link[nontarget]{pattern.search}}
	\code{\link[nontarget]{combine}}
}

