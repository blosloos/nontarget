\name{nontarget-package}
\alias{nontarget-package}
\alias{nontarget}
\docType{package}
\title{
Detecting Isotope, Adduct and Homologue Relations in LC-MS Data.
}
\description{
Grouping of peaks in a HRMS data set for (1) isotopic pattern relations and (2) different adducts of the same molecule; detection of (3) homologue series. 
Isotopic pattern and adduct groups can then be related to their (unknown) candidate chemical component, with homologue series information attached.
Includes various plotting and filtering functions for e.g. mass defects, frequent m/z distances, components vs. non-components, adduct frequencies.  
}
\details{
\tabular{ll}{
Package: \tab nontarget\cr
Type: \tab Package\cr
Version: \tab 1.8\cr
Date: \tab 2015-02-09\cr
License: \tab GPL-3\cr
}
Screens a HRMS data set for peaks related by (1) isotopic patterns and/or (2) different adducts of the same molecule and/or (3) being part of 
a homologue series, including various plausibility checks. The resulting isotopic pattern groups and adduct groups can then be combined to components,
with each component tagged if being part of a homologue series. This does not require prior knowledge about the chemical nature of the components assigned.


Includes various plotting functions, such as (a) m/z vs. RT vs. mass defect, (b) mass defect vs. detected isotope m/z increments, (c) adduct frequencies and their
intensity distributions, (d) relations among peaks within single isotope/adduct groups and within single components and (e) homologue series within RT vs. m/z plots.
Allows filtering HRMS data for mass defects, satellite peaks, frequent m/z distances and components vs. non-components.
Lists of most-common adducts and isotopes are provided or may be user-defined.


Requires HRMS centroid peak lists as input, i.e., a dataframe or matrix with values of (a) m/z, (b) intensity and (c) retention time (RT) per peak.
In addition, tolerances for m/z, RT and uncertainties in peak intensity must be defined by the user.
}
\author{
Martin Loos

Maintainer: Martin Loos <Martin.Loos@eawag.ch>

}
\references{
Loos, M., Hollender, J., Schymanski, E., Ruff, M., Singer, H., 2012.
Bottom-up peak grouping for unknown identification from high-resolution mass spectrometry data.
ASMS 2012 annual conference Vancouver, oral session Informatics: Identification.
}

\keyword{ package }

\seealso{
	
	Detecting isotope pattern groups:
	\code{\link[nontarget]{peaklist}}
	\code{\link[nontarget]{make.isos}}
	\code{\link[nontarget]{pattern.search}}
	\code{\link[nontarget]{pattern.search2}}
	\code{\link[nontarget]{plotisotopes}}	
	\code{\link[nontarget]{plotdefect}}
	\code{\link[enviPat]{isotopes}}
	\code{\link[enviPat]{resolution_list}}
	
	Detecting adduct groups:
	\code{\link[nontarget]{peaklist}}
	\code{\link[nontarget]{adduct.search}}	
	\code{\link[nontarget]{plotadduct}}	
	\code{\link[enviPat]{adducts}}
	
	Detecting homologue series:
	\code{\link[nontarget]{peaklist}}	
	\code{\link[nontarget]{homol.search}}
	\code{\link[nontarget]{plothomol}}

	
	On combining groups to components:
	\code{\link[nontarget]{combine}}	
	\code{\link[nontarget]{plotisotopes}}	
	\code{\link[nontarget]{plotcomp}}
	\code{\link[nontarget]{ms.filter}}
	
	On filtering and plotting:
	\code{\link[nontarget]{rm.sat}}
	\code{\link[nontarget]{plotall}}	
	\code{\link[nontarget]{plotgroup}}	
	\code{\link[nontarget]{ms.filter}}
	\code{\link[nontarget]{plotdiff}}
	\code{\link[nontarget]{deter.iso}}  	
	

  }


\examples{
\donttest{
######################################################
# (0) load required data: ############################
# (0.1) HRMS peak list & remove satelite peaks: ######
data(peaklist);
peaklist<-rm.sat(peaklist,dmz=0.3,drt=0.1,intrat=0.015,spar=0.8,corcut=-1000,plotit=TRUE);
peaklist<-peaklist[peaklist[,4],1:3];
# (0.2) list of adducts - package enviPat ############
data(adducts);
# (0.3) list of isotopes - package enviPat ###########
data(isotopes);
######################################################
# (1) run isotope pattern grouping ###################
# (1.1) define isotopes and charge argument ##########
iso<-make.isos(isotopes,
	use_isotopes=c("13C","15N","34S","37Cl","81Br","41K","13C","15N","34S","37Cl","81Br","41K"),
	use_charges=c(1,1,1,1,1,1,2,2,2,2,2,2))
# (1.2) run isotope grouping #########################
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
# (1.3) plot results #################################
plotisotopes(pattern);
plotdefect(pattern,elements=c("N"));
######################################################
# (2.1) run grouping of peaks for different adducts ##
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
# (2.2) plot results #################################
plotadduct(adduct);
######################################################
# (3) show single pattern group and its relation #####
# to adduct groups ###################################
plotall(pattern,adduct);
plotgroup(pattern,adduct,groupID=1,massrange=10,allmass=FALSE);
######################################################
# (4.1) Screen for homologue series ##################
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
# (4.2) Plot results #################################
plothomol(homol,xlim=FALSE,ylim=FALSE,plotlegend=TRUE);
######################################################
# (5.1) Combine grouping results to components #######
comp<-combine(
	pattern,
	adduct,
	homol,
	dont=FALSE,
	rules=c(TRUE,FALSE,FALSE)
);
# (5.2) plot results ################################# 
plotisotopes(comp);
plotcomp(comp,compoID=1,peakID=FALSE);
######################################################
# (6) Select data from interactive plot ##############
# ms.filter( component=comp,x="mz",y="dm",xlim=FALSE,
# ylim=FALSE,rm.comp=TRUE,plot.comp=TRUE,rm.noncomp=FALSE,
# select.polygon="inside",res=100,filter.for="raw" );
######################################################

}
}





















