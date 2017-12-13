\name{adduct.search}
\alias{adduct.search}
\title{
Detecting and grouping adduct m/z relations among peaks in a HRMS dataset
}
\description{
Algorithm for detecting m/z differences among peaks that may correspond to m/z differences among different adducts.
}
\usage{
adduct.search(peaklist, adducts, rttol = 0, mztol = 2,ppm = TRUE, 
use_adducts = c("M+H", "M+K", "M+Na"), ion_mode = "positive", get_pairs = FALSE)
}


\arguments{
  \item{peaklist}{
	Dataframe of HRMS peaks with three numeric columns for (a) m/z, (b) intensity and (c) retention time, such as \code{\link[nontarget]{peaklist}}.
}
  \item{adducts}{
	Data.frame \code{\link[enviPat]{adducts}} or equivalent.
}
  \item{rttol}{
	Retention time tolerance. Units as given in column 3 of peaklist argument, e.g. [min].
}
  \item{mztol}{
	m/z tolerance setting: +/- value by which the m/z of a measured peak may vary from its expected m/z value. If parameter \code{ppm=TRUE} (see below) given in ppm, otherwise, if \code{ppm=FALSE}, in absolute m/z [u]. 
}
  \item{ppm}{
	Should \code{mztol} be set in ppm (TRUE) or in absolute m/z (FALSE)
}
  \item{use_adducts}{
	Vector of adducts to be screened for. Corresponds to names in the first column of \code{\link[enviPat]{adducts}}, thus referring to equations from the second column of \code{\link[enviPat]{adducts}}
	to be used for calculating adduct m/z differences.
}
  \item{ion_mode}{
	\code{"positive"} or \code{"negative"}.
}
  \item{get_pairs}{
	enviMass output, please ignore.
}
}
\details{
Given a peak from the peaklist, the \code{adduct.search} algorithm screens within tolerances \code{mztol} and \code{rttol} whether any other peaks may correspond to this one peak via adduct m/z differences.
More precisely, the one peak m/z is reset to all possible candidate molecular mass values (M; uncharged, non-adduct). The latter are then used to calculate for all other candidate adduct peaks,
which, if found, are subsequently grouped.

For example, consider \code{use_adducts=c("M+H", "M+K")}. Given the m/z-value of the one peak, two other peaks with 
((m/z*z("M+H")-X("M+H"))/z("M+K"))+X("M+K") and ((m/z*z("M+K")-X("M+K"))/z("M+H"))+X("M+H") are searched for. The peak found for
the first term (i.e. with "M+H" being the candidate adduct of the one peak) leads to one group of associated adduct peaks (M+H<->M+K). Another adduct
peak (i.e. with "M+K" being the candidate adduct of the one peak) would lead to a second group of associated adduct peaks (M+K<->M+H). 
Logically, larger adduct groups than the one exemplified can be present, if argument "use_adducts" allows for it (e.g. M+H<->M+K,M+H<->M+Na,M+Na<->M+K).

For clarification, \code{mztol} states the maximum m/z deviation of a measured peak from its true value, i.e., the theoretical mass-to-charge ratio of the (often unknown) analyte adduct measured. The latter true value thus ranges +/-\code{mztol} from
the measured value, leading to a lower and an upper m/z bound for this true value. These bounds are then modified by pairwise adduct m/z differences, leading to new bounds for the true value at other m/z positions. In turn, the values of measured peaks at 
exactly these positions can again deviate by +/-\code{mztol} from these bounds, which are hence adapted accordingly at these positions for a final search window. 
This entails: (a) the bounds are calculated from measured instead of true values (using \code{ppm==TRUE}, bound differences are assumed negligible between utilizing measured or true values) and 
(b) the final search window is larger than 2*\code{mztol} and can therefore lead to contradictory assignments (i.e., peaks B and C can be assigned to peak A for different adduct m/z differences within \code{mztol} but peaks B and C cannot be paired 
within \code{mztol}).
}
\value{

List of type adduct with 5 entries 

\item{adduct[[1]]}{\code{Adducts}. Dataframe with peaks (\code{mass},\code{intensity},\code{rt},\code{peak ID}) and their adduct
relations (\code{to ID},\code{adduct(s)},\code{mass tolerance},\code{charge level}) within adduct groups (\code{group ID},\code{interaction level}).}
\item{adduct[[2]]}{\code{Parameters}. Parameters used.}
\item{adduct[[3]]}{\code{Peaks in adduct groups}.Dataframe listing all peaks (\code{peak IDs}) for an adduct group 
(\code{group ID}) and the individual adducts found in that group (\code{adducts}).}
\item{adduct[[4]]}{\code{Number of adducts}. Counts of hits per adduct over all adduct groups found.}
\item{adduct[[5]]}{\code{Overlaps}. Count on how many peaks were assigned to be two different adducts}

}
\author{
Martin Loos
}
\note{
Peak IDs refer to the order in which peaks are provided. 
Different IDs exist for adduct groups, isotope pattern groups, grouped homologue series (HS) peaks 
and homologue series cluster. Yet other IDs exist for the individual components (see note section of \code{\link[nontarget]{combine}}).

The same peak may appear as different adducts in column \code{adduct[[1]][,7]}, indicating a conflict in assigning the correct adduct.
Beware, some adduct combinations from \code{\link[enviPat]{adducts}} may lead to the same results (e.g. M+H<->M+Na vs M+3H<->M+3Na).

}

\seealso{
	\code{\link[nontarget]{rm.sat}}
	\code{\link[enviPat]{adducts}}
	\code{\link[nontarget]{peaklist}}
	\code{\link[nontarget]{plotadduct}}	
	\code{\link[nontarget]{combine}}
	\code{\link[nontarget]{plotgroup}}
}
\examples{
\donttest{
######################################################
# load required data: ################################
# HRMS peak list: ####################################
data(peaklist)
# list of adducts ####################################
data(adducts) 
######################################################
# run grouping of peaks for different adducts ########
# of the same candidate molecule #####################
adduct<-adduct.search(
  peaklist,
  adducts,
  rttol=0.05,
  mztol=3,
  ppm=TRUE,
  use_adducts=c("M+K","M+H","M+Na","M+NH4"),
  ion_mode="positive"
);
# plot results #######################################
plotadduct(adduct);
######################################################
}
}
