\name{pattern.search2}
\alias{pattern.search2}

\title{
Detecting and grouping isotope m/z relations among LC-HRMS centroid peaks, based on quantized reference data
}
\description{
Algorithm for grouping isotope pattern centroids of chemical components by querying quantized simulation data
}
\usage{
pattern.search2(peaklist,quantiz,mztol=2,ppm=TRUE,inttol=0.5,rttol=0.3,
use_isotopes=c("13C","37Cl","15N","81Br","34S","18O"),use_charges=c(1,2),
use_marker=TRUE,quick=FALSE,isotopes, get_pairs = FALSE, get_matches = FALSE)
}


\arguments{
  \item{peaklist}{
	Dataframe of HRMS peaks with three numeric columns for (a) m/z, (b) intensity and (c) retention time, such as \code{\link[nontarget]{peaklist}}.
}
  \item{quantiz}{
	Quantized simulation data of feasible centroid-centroid relations as provided by package \code{nontargetData}.
}
  \item{mztol}{
	m/z tolerance setting: value by which the m/z of a peak may vary from its expected value. If parameter \code{ppm=TRUE} (see below) given in ppm, otherwise, if \code{ppm=FALSE}, in absolute m/z [u].
}
  \item{ppm}{
	Should \code{mztol} be set in ppm (\code{TRUE}) or in absolute m/z (\code{FALSE}).
}
  \item{inttol}{
	Intensity tolerance setting = fraction by which peak intensities may vary; e.g., if set to 0.2, a peak with expected intensity 10000 may range in between 8000 and 12000.	
}
  \item{rttol}{
	+/- retention time tolerance. Units as given in column 3 of \code{peaklist} argument, e.g. [min].
}
  \item{use_isotopes}{
	Restrict query to certain isotopes dominating centroid relations; set to \code{FALSE} to use all available isotopes.
}
  \item{use_charges}{
	Vector of signed integers. Restrict query to certain charges z; set to \code{FALSE} to use all charge states. 
}  
  \item{use_marker}{  
	Query for marker peaks, \code{FALSE} or \code{TRUE}?
}  
  \item{quick}{  
	Continue if query finds first hit? Speeds up, but leaves resulting information on underlying isotopes incomplete.  
} 
  \item{isotopes}{ 
	Dataframe of relevant isotopes as provided by package \code{enviPat}; used for checking user inputs.
}  
  \item{get_pairs}{
	enviMass output, please ignore.
}
  \item{get_matches}{
	enviMass output, please ignore.
}
}

\details{
As alternative to rule-based \code{\link[nontarget]{pattern.search}}, differences among measured centroids (\code{peaklist}) are queried to match those of compressed (=quantized) simulation data within bounds 
of measurement tolerances and the quantization distortion. Hence, in comparion to \code{\link[nontarget]{pattern.search}}, this approach accounts for centroid mass shifts induced by peak profile 
interferences prevalent at even high m/z resolution.

To derive the quantized data, isotope pattern centroids of several million organic molecular formulas from the PubChem database were calculated for various classes of adducts. 
Molecular formulas were filtered to be unique and only to contain C, H, O, N, Cl, Br, K, Na, S, Si, F, P and/or I. 
The resulting >250 million centroid pairs from individual patterns were then categorized for their dominant isotopologues, charge and
the possible presence of another centroid of higher intensity than that of the pair (=marker peak).
Within these categories, data on centroid pair (a) m/z, (b) m/z differences, (c) intensity ratios and (d) marker m/z was quantized by a 
recursive partitioning procedure.
The resulting compressed data representation was extended by nearest neigbour estimates in the above dimensions (a) to (d) to 
account for queries with molecular formulas possibly not present in the PubChem set. 
Internally, the quantized simulation data is queried by a tree-like space-partitioning structure for hyperrectangles, while centroids from \code{peaklist} are restructured into kd-trees. 
}

\note{
Peak IDs refer to the order in which peaks are provided. 


If you do not find quantized simulation data for your instrument in package \code{nontargetData} and you can provide resolution=f(m/z) information: contact maintainer.
}

\section{Warning}{
Acceptable outcomes strongly depend on appropriate parametrization of the algorithm and using the correct \code{quantiz} data set from package \code{nontargetData}. 
Using overly large values for \code{rttol} and/or \code{mztol} may lead to slow execution. 
}

\value{
List of type pattern with 12 entries

\item{pattern[[1]]}{\code{Patterns}. Dataframe with peaks (\code{mass},\code{intensity},\code{rt},\code{peak ID}) and their 
isotope pattern relations (\code{to ID},\code{isotope(s)},\code{mass tolerance} (deprecated),\code{charge level}) within
isotope pattern groups (\code{group ID},\code{interaction level} (deprecated)).}
\item{pattern[[2]]}{\code{Parameters}. Parameters used.}
\item{pattern[[3]]}{\code{Peaks in pattern groups}. Dataframe listing all peaks (\code{peak IDs}) per isotope pattern group (\code{group ID}) at the given z-level(s) (\code{charge level}).}
\item{pattern[[4]]}{\code{Atom counts}. Deprecated.}
\item{pattern[[5]]}{\code{Count of pattern groups}. Number of isotope pattern groups found on the different z-levels used.}
\item{pattern[[6]]}{\code{Removals by rules}. Deprecated.}
\item{pattern[[7]]}{\code{Number of peaks with pattern group overlapping}. Deprecated}
\item{pattern[[8]]}{\code{Number of peaks per within-group interaction levels}.}
\item{pattern[[9]]}{\code{Counts of isotopes}. Number of times a m/z isotope difference was detected (raw measure / number of isotope pattern groups)}
\item{pattern[[10]]}{\code{Elements}. Elements used via argument iso derived by \code{\link[nontarget]{make.isos}}.}
\item{pattern[[11]]}{\code{Charges}. z-levels used.}
\item{pattern[[12]]}{\code{Rule settings}. Deprecated.}

}
\author{
Martin Loos
}
\seealso{
	\code{\link[nontarget]{rm.sat}}
	\code{\link[nontarget]{peaklist}}
	\code{\link[nontarget]{plotisotopes}}	
	\code{\link[nontarget]{plotdefect}}
	\code{\link[nontarget]{combine}}
	\code{\link[nontarget]{plotgroup}}
	\code{\link[nontarget]{pattern.search}}
}

\examples{
\donttest{
######################################################
# load HRMS centroid list: ###########################
data(peaklist)
# load isotope data ##################################
data(isotopes)
# load quantized simulation data #####################
data(OrbitrapXL_VelosPro_R60000at400_q)
######################################################
# run isotope pattern grouping #######################
# save the list returned as "pattern" ################
pattern<-pattern.search2(
	peaklist,
	OrbitrapXL_VelosPro_R60000at400_q,
	mztol=2, 
	ppm=TRUE,
	inttol=0.5,
	rttol=0.3,
	use_isotopes=FALSE,
	use_charges=FALSE,
	use_marker=TRUE,
	quick=FALSE,
	isotopes
)
names(pattern);
######################################################
}
}
