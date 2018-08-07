\name{BCG}
\alias{BCG}
\docType{data}
\title{Dataset on the Effectiveness of the BCG Vaccine for Preventing Tuberculosis
}
\description{
  This dataset includes 13 studies on the effectiveness of the
  Bacillus Calmette-Guerin (BCG) vaccine for preventing tuberculosis
  (see van Houwelingen, Arends, & Stijnen (2002) for details).
}
\usage{data(BCG)}

\details{
A list of data with the following structure:
\describe{
  \item{Trial}{Number of the trials}
  \item{Author}{Authors of the original studies}
  \item{Year}{Year of publication}
  \item{VD}{Vaccinated group with disease}
  \item{VWD}{Vaccinated group without disease}
  \item{NVD}{Not vaccinated group with disease}
  \item{NVWD}{Not vaccinated group without disease}
  \item{Latitude}{Geographic latitude of the place where the study was done}
  \item{Allocation}{Method of treatment allocation}
  \item{ln_OR}{Natural logarithm of the odds ratio: log((VD/VWD)/(NVD/NVWD))} 
  \item{v_ln_OR}{Sampling variance of ln_OR: 1/VD+1/VWD+1/NVD+1/NVWD}
  \item{ln_Odd_V}{Natural logarithm of the odds of the vaccinated group:
	log(VD/VWD)}
  \item{ln_Odd_NV}{Natural logarithm of the odds of the not vaccinated group: log(NVD/NVWD)}  
  \item{v_ln_Odd_V}{Sampling variance of ln_Odd_V: 1/VD+1/VWD}
  \item{cov_V_NV}{Sampling covariance between ln_Odd_V and ln_Odd_NV: It
  is always 0}
  \item{v_ln_Odd_NV}{Sampling variance of ln_Odd_NV: 1/NVD+1/NVWD}  
}
}

\source{
Colditz, G. A., Brewer, T. F., Berkey, C. S., Wilson, M. E., Burdick, E., Fineberg, H. V., & Mosteller, F. (1994). Efficacy of BCG vaccine in the prevention of tuberculosis: Meta-analysis of the published literature. \emph{Journal of the American Medical Association}, \bold{271}, 698--702.
}
\references{
Berkey, C. S., Hoaglin, D. C., Mosteller, F., & Colditz, G. A. (1995). A random-effects regression model for meta-analysis. \emph{Statistics in Medicine}, \bold{14}, 395--411.

van Houwelingen, H. C., Arends, L. R., & Stijnen, T. (2002). Advanced methods in meta-analysis: Multivariate approach and meta-regression. \emph{Statistics in Medicine}, \bold{21}, 589--624.

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}

\examples{
data(BCG)

## Univariate meta-analysis on the log of the odds ratio
summary( meta(y=ln_OR, v=v_ln_OR, data=BCG,
              x=cbind(scale(Latitude,scale=FALSE),
              scale(Year,scale=FALSE))) )

## Multivariate meta-analysis on the log of the odds
## The conditional sampling covariance is 0
bcg <- meta(y=cbind(ln_Odd_V, ln_Odd_NV), data=BCG,
            v=cbind(v_ln_Odd_V, cov_V_NV, v_ln_Odd_NV))
summary(bcg)

plot(bcg)
}
\keyword{datasets}
