\name{Nam03}
\alias{Nam03}
\docType{data}
\title{Dataset on the Environmental Tobacco Smoke (ETS) on children's health
}
\description{
  This dataset includes 59 studies reported by Nam, Mengersen, and
  Garthwaite (2003) on the potential health effects among children exposed to
environmental tobacco smoke (ETS), or passive smoking. The effect sizes
are the log odds ratios of asthma and lower respiratory disease (LRD).
}
\usage{data(Nam03)}

\details{
A list of data with the following structure:
\describe{
  \item{ID}{Study identification number.}
  \item{Size}{Total number of valid subjects in the study.}
  \item{Age}{Mean age of participants.}
  \item{Year}{Year of publication.}
  \item{Country}{Country code.}
  \item{Smoke}{Source of ETS.}
  \item{Adj}{Whether the reported odds ratio is adjusted for covariates.}
  \item{Asthma_logOR}{Log odds ratio of asthma.}
  \item{LRD_logOR}{Log odds ratio of lower respiratory disease.}
  \item{Asthma_v}{Sampling variance of Asthma_logOR.}
  \item{AsthmaLRD_cov_05}{Sampling covariance between Asthma_logOR and
	LRD_logOR by assuming a correlation of 0.5}
  \item{LRD_v}{Sampling variance of LRD_logOR.}  
}
}

\source{
  Nam, I.-S., Mengersen, K., & Garthwaite, P. (2003). Multivariate meta-analysis. \emph{Statistics in Medicine}, \bold{22}(14), 2309-2333. https://doi.org/10.1002/sim.1410
}

\examples{
data(Nam03)
}
\keyword{datasets}
