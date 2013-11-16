\name{Jaramillo05}
\alias{Jaramillo05}
\docType{data}
\title{Dataset from Jaramillo, Mulki & Marshall (2005)
}
\description{A dataset of the relationship between organizational
  commitment (OC) and salesperson job performance (JP) from Jaramillo, Mulki & Marshall (2005).
}
\usage{data(Jaramillo05)}
\format{
  A data frame with 61 observations on the following 10 variables.
  \describe{
    \item{\code{Author}}{a character vector of study}
    \item{\code{Sample.size}}{sample size of the study}
    \item{\code{Sales}}{sample type; either "mixed", "nonsales" or "sales"}
    \item{\code{Country}}{a character vector of country of study}
    \item{\code{IDV}}{Hofstede's (1997) individualism index}
    \item{\code{OC.scale}}{scale of OC; either "Porter or Mowday",
	  "Meyer" or "other"}
    \item{\code{OC.alpha}}{Coefficient alpha of organizational commitment}
    \item{\code{JP.alpha}}{Coefficient alpha of job performance}
    \item{\code{r}}{correlation between organizational commitment and
	  job performance}
    \item{\code{r.v}}{sampling variance of r}
  }
}

\source{Jaramillo, F., Mulki, J. P., & Marshall, G. W. (2005). A meta-analysis of the relationship between organizational commitment and salesperson job performance: 25 years of research. \emph{Journal of Business Research}, \bold{58(6)}, 705-714. doi:10.1016/j.jbusres.2003.10.004
}

\examples{
## Research question 4.4.1
summary(meta(r, r.v, data=Jaramillo05))

## Research question 4.4.2
## Select cases with either "sales" or "nonsales"
Sales.df <- subset(Jaramillo05, Sales \%in\% c("sales", "nonsales"))

## Create a predictor with 1 and 0 when they are "sales" or "nonsales", respectively
predictor <- ifelse(Jaramillo05$Sales=="sales", yes=1, no=0)

## Mixed-effects meta-analysis
summary( meta(y = r, v = r.v, x = predictor, data = Jaramillo05) )

## Research question 4.4.3
summary(meta(r, r.v, x=IDV, data=Jaramillo05))
}
\keyword{datasets}
