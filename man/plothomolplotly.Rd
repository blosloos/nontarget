\name{plothomolplotly}
\alias{plothomolplotly}
\title{
  Creates interactive visualization from the output of homol.search() of retention time (RT) vs. m/z
}
\description{
  Given results from \code{\link[nontarget]{homol.search}}, an interactive plot of homologue series (HS) is generated.
  The homologue series are color coded in groups of common m/z increments. Hovering over the HS reveals further information of the peaks.
}
\usage{
  plothomolplotly(homol)
}

\arguments{
  \item{homol}{
    List of type homol produed by \code{\link[nontarget]{homol.search}}.
  }
  
}
\author{
  Martin Loos; Diana Masch
}

\seealso{
  \code{\link[nontarget]{homol.search}}
}
\examples{
  \donttest{
    data(peaklist);
    data(isotopes)
    homol<-homol.search(
      peaklist,
      isotopes,	
      elements=c("C","H","O"),
      use_C=TRUE,
      minmz=5,
      maxmz=120,
      minrt=2,
      maxrt=2,
      ppm=TRUE,
      mztol=3.5,
      rttol=0.5,
      minlength=3,
      mzfilter=FALSE,
      vec_size=3E6,
      spar=.45,
      R2=.98,
      plotit=FALSE
    )
    plothomolplotly(homol);
  }
}

