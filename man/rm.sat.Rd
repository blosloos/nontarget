\name{rm.sat}
\alias{rm.sat}

\title{
Removal of satellite peaks from FT-MS peak lists
}
\description{
Brute force method to remove satellite peak from a FT-HRMS peak list.
}
\usage{
rm.sat(peaklist, dmz = 0.3, drt = 0.3, intrat = 0.01, spar = 0.8, 
corcut = 0.8, plotit = TRUE)
}

\arguments{
  \item{peaklist}{
	Dataframe of HRMS peaks with three numeric columns for (a) m/z, (b) intensity and (c) retention time, such as \code{\link[nontarget]{peaklist}}.
	}
  \item{dmz}{
	m/z window around a parent peak within which satellite peaks are searched for.
	}
  \item{drt}{
	Retention time window around a parent peak within which satellite peaks are searched for.
	}
  \item{intrat}{
	Intensity ratio between satellite peak/associated parent peak below which the former are removed. 
	}
  \item{spar}{
	\code{spar} argument used in R-function smooth spline. See details and \code{?smooth.spline}.
	}
  \item{corcut}{
	Correlation coefficient above which symmetrical peaks are marked as satellite peaks. To disable, set to \code{-1000} See details. 
	}
  \item{plotit}{
	Plot results?
}
}
\details{
"Parent" peak refers to a peak having associated satellite peaks as artifacts from FT calculations.

\code{\link[nontarget]{rm.sat}} screens, along decreasing intensity, peaks for having other peaks within ranges set by arguments \code{dmz}, \code{drt} and \code{intrat}.
If present, the latter are marked as satellite peaks and are subsequently excluded from further screening within \code{\link[nontarget]{rm.sat}}.

In addition, arguments \code{spar} and \code{corcut} evaluate the symmetry of satellite peaks around the parent peak (i.e. below and above the parent peak m/z), if
enough peaks around a parent peak within ranges set by arguments \code{dmz}, \code{drt} and \code{intrat} are found (here: at least 8 peaks, 4 above and 4 below the parent peak m/z). 
Two splines are fitted by R function \code{smooth.spline}, one to those peaks above and one to those peaks below the parent peak m/z.
If the splines are symmetric (i.e. correlated with each other, see argument corcut), the associated peaks are termed satellites. This approach has not yet faced 
validation and is highly dependent on the peak-picking algorithm.
}
\value{
A dataframe with four columns. The first three columns are identical to those of argument \code{peaklist}.
The fourth columns marks potential satellite peaks with \code{FALSE}, the other peaks with \code{TRUE} (see example).
}
\author{
Martin Loos
}
\note{
Not removing satellite peaks may lead to undesirable artifacts when screening for isotope pattern and adduct relations using \code{\link[nontarget]{pattern.search}} 
and \code{\link[nontarget]{adduct.search}}, respectively. For example, consider a satellite peaks having a slightly larger m/z than its monoisotopic parent peaks. 
Then, a m/z difference from a 13C isotope between monoisotopic parent and M+1 peak often leads to a 15N isotope difference between satellite and M+1 peak. This artifact 
causes bogus isotope pattern groups (with the satellite peak assigned the monoisotopic peak in this example), group overlaps (see \code{\link[nontarget]{pattern.search}})
and interfering peaks in components (see \code{\link[nontarget]{combine}}).

Still, given the brute approach of \code{\link[nontarget]{rm.sat}}, there is no guarantee that all peaks removed are indeed satellite peaks. As an alternative, one may
filter for peaks with an overly short eluation time/scan number or use data from MS devices that are less prone to produce satellite peaks.
}

\seealso{
	\code{\link[nontarget]{peaklist}}
}
\examples{
\donttest{
data(peaklist);
peaklist<-rm.sat(peaklist,dmz=0.3,drt=0.1,intrat=0.015,spar=0.8,corcut=-1000,plotit=TRUE);
peaklist<-peaklist[peaklist[,4],1:3];
}
}
