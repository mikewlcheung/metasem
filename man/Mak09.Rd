\name{Mak09}
\alias{Mak09}
\docType{data}
\title{Eight studies from Mak et al. (2009)}
\description{
  Eight studies from Mak et al. (2009) were reported by Cheung et al. (2012).
}

\usage{data(Mak09)}
\format{
  A data frame with 8 observations on the following 10 variables.
  \describe{
    \item{\code{Study}}{a character vector of study}
    \item{\code{type}}{a character vector}
    \item{\code{AF.BP}}{a numeric vector}
    \item{\code{Tot.BP}}{a numeric vector}
    \item{\code{AF.non.BP}}{a numeric vector}
    \item{\code{Tot.non.BP}}{a numeric vector}
    \item{\code{yi}}{a numeric vector}
    \item{\code{vi}}{a numeric vector}
    \item{\code{age.mean}}{a numeric vector}
    \item{\code{study.duration}}{a numeric vector}
  }
}

\source{
  Mak, A., Cheung, M. W.-L., Ho, R. C. M., Cheak, A. A. C., & Lau,
  C. S. (2009). Bisphosphonate and atrial fibrillation: Bayesian
  meta-analyses of randomized controlled trials and observational
  studies. \emph{BMC Musculoskeletal Disorders},
  \bold{10(113)}. doi:10.1186/1471-2474-10-113 Available at \url{https://bmcmusculoskeletdisord.biomedcentral.com/articles/10.1186/1471-2474-10-113}.
}
\references{
 Cheung, M. W.-L., Ho, R. C. M., Lim, Y., & Mak, A. (2012). Conducting a
 meta-analysis: Basics and good practices. \emph{International Journal
   of Rheumatic Diseases}, \bold{15(2)}, 129-135. doi: 10.1111/j.1756-185X.2012.01712.x
}
\examples{
## Random-effects meta-analysis
( meta1 <- summary(meta(y=yi, v=vi, data=Mak09, I2=c("I2q", "I2hm"))) )

## Convert the estimates back into odds ratio 
OR <- with(coef(meta1), exp(c(Estimate[1], lbound[1], ubound[1])))
names(OR) <- c("Estimate in OR", "lbound in OR", "ubound in OR")
OR
  
## Mixed-effects meta-analysis with mean age as a predictor
summary( meta(y=yi, v=vi, x=age.mean, data=Mak09) )
}
\keyword{datasets}
