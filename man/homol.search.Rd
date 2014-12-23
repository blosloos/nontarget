\name{homol.search}
\alias{homol.search}

\title{
Detecting homologue series in HRMS data.
}
\description{
An algorithm for unsupervised detection of homologue series in HRMS data.
}
\usage{
homol.search(peaklist,isotopes,elements=c("C","H","O"),charge=c(1,2),use_C=TRUE,
minmz=5,maxmz=60,minrt=0.5,maxrt=2,ppm=TRUE,mztol=3.5,rttol=0.5,minlength=4,
mzfilter=FALSE,vec_size=1E6)
}


\arguments{
  \item{peaklist}{
	Dataframe of HRMS peaks with three numeric columns for (a) m/z, (b) intensity and (c) retention time, such as \code{\link[nontarget]{peaklist}}.
}
  \item{isotopes}{
	Dataframe \code{\link[enviPat]{isotopes}}
}

  \item{elements}{ 
  Chemical elements considered for the deviation within a homologue series, e.g. c("C","H") for alkane chains. Used to restrict search.
}

  \item{charge}{
	Possible MS ionization charge of a homologue series. Used to restrict search.
}

  \item{use_C}{
	For \code{elements}: take element ratio to C-atoms into account? Used to restrict search.
}


  \item{minmz}{
	Defines the lower limit of the m/z window to look for homologue series peaks in, relative to the m/z of the one peak to search from. Given in absolute m/z values [u]. 
}

  \item{maxmz}{
	Defines the upper limit of the m/z window to look for homologue series peaks in, relative to the m/z of the one peak to search from. Given in absolute m/z values [u]. 
}

  \item{minrt}{
	Defines the lower limit of the retention time (RT) window to look for other homologue peaks, relative to the RT of the one peak to search from. 
}

  \item{maxrt}{
	Defines the upper limit of the retention time (RT) window to look for other homologue peaks, relative to the RT of the one peak to search from. 
}

  \item{ppm}{
	Should \code{mztol} be set in ppm (\code{TRUE}) or in absolute m/z [u] (\code{FALSE})?
}

  \item{mztol}{
	m/z tolerance setting: value by which the m/z of a peak may vary from its expected value. If parameter \code{ppm=TRUE} (see below) given in ppm, otherwise, if \code{ppm=FALSE}, in absolute m/z [u]. 
}

  \item{rttol}{
	Retention time (RT) tolerance by which the RT between two adjacent pairs of a homologue series is allowed to differ. Units as given in column 3 of peaklist argument, e.g. [min].
}

  \item{minlength}{
	Minimum number of peaks in a homologue series.
}

  \item{mzfilter}{
	Vector of numerics. Can be used to filter output for homologue series with certain differences in m/z along the series, given the tolerances in \code{mztol}.
}

  \item{vec_size}{
	Vector size, used for debugging. Ignore unless the relevant error message is printed.
}

}

\note{
Input \code{isotopess} is needed to limit intermediate numbers of m/z differences to screen over, based on changes in mass defect. 
Similarly, intermediate numbers are limited by the retention time and m/z windows defined by the \code{minrt / maxrt} and  \code{minmz / maxmz} arguments. The latter are always set relative to
the individual RT and m/z values of the peaks to be searched from. 

Peak IDs refer to the order in which peaks are provided. Different group IDs exist for adduct groups, isotope pattern groups and homologue series peak 
groups. Yet other IDs exist for the individual components (see note section of \code{\link[nontarget]{combine}}).

Here, IDs of homologue series group are given both in the function output \code{homol[[1]]} and \code{homol[[3]]}, with one homologue series stating one group of interrelated peaks.
}

\section{Warning}{
The \code{rttol} argument of \code{\link[nontarget]{homol.search}} must not be mixed with that of \code{\link[nontarget]{pattern.search}} or \code{\link[nontarget]{pattern.search2}}.
}
\value{

List of type homol with 5 entries 

\item{homol[[1]]}{\code{Homologue Series}. Dataframe with peaks (\code{mass},\code{intensity},\code{rt},\code{peak ID}) and their homologue series relations (\code{to ID},\code{m/z increment},\code{RT increment}) within
groups of homologue series (\code{group ID},\code{series level}).}
\item{homol[[2]]}{\code{Parameters}. Parameters used.}
\item{homol[[3]]}{\code{Peaks in homologue series}. Dataframe listing all peaks (\code{peak IDs}) per homologue series (\code{group ID}), the underlying mean m/z & RT increments 
(\code{m/z increments}, \code{RT increments}) 
and the minimum and maximum RT changes between individual peaks of the series.}
\item{homol[[4]]}{\code{m/z restrictions used}. See function argument \code{mzfilter}.}
\item{homol[[5]]}{\code{Peaks per level}. List of peak IDs per level in the individual series.}
\item{homol[[6]]}{\code{Single relations}. Peak to peak relations of series.}
}

\author{
Martin Loos
}
\seealso{
	\code{\link[nontarget]{rm.sat}}
	\code{\link[enviPat]{isotopes}}
	\code{\link[nontarget]{peaklist}}
	\code{\link[nontarget]{plothomol}}
}


\examples{
\donttest{
data(peaklist);
data(isotopes)
homol<-homol.search(
		peaklist,
		isotopes,
		elements=c("C","H","O"),
		charge=c(1,2),
		use_C=TRUE,
		minmz=5,
		maxmz=60,
		minrt=0.5,
		maxrt=2,
		ppm=TRUE,
		mztol=3.5,
        rttol=0.5,
		minlength=4,
		mzfilter=FALSE,
		vec_size=1E6
)
plothomol(homol);
}
}