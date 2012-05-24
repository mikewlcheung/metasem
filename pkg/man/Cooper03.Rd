\name{Cooper03}
\alias{Cooper03}
\docType{data}
\title{Selected effect sizes from Cooper et al. (2003)}
\description{
  Fifty-six studies from Cooper et al. (2003) were reported by Konstantopoulos (2011).
}
\usage{data(Cooper03)}
\format{
}
\details{
    The variables are:
  \describe{
  \item{District}{District ID}
  \item{Study}{Study ID}
  \item{y}{Effect size}
  \item{v}{Sampling variance}
  \item{Year}{Year of publication}
}

}
\source{
  Cooper, H., Valentine, J. C., Charlton, K., & Melson, A. (2003). The Effects of Modified School Calendars on Student Achievement and on School and Community Attitudes. \emph{Review of Educational Research}, \bold{73(1)}, 1-52. doi:10.3102/00346543073001001
}
\references{
  Konstantopoulos, S. (2011). Fixed effects and variance components
  estimation in three-level meta-analysis. \emph{Research Synthesis
  Methods}, \bold{2}, 61-76. doi:10.1002/jrsm.35
}
\examples{
data(Cooper03)

## No predictor
summary( meta3(y=y, v=v, cluster=District, data=Cooper03) )

## Show all heterogeneity indices and their 95\% confidence intervals
summary( meta3(y=y, v=v, cluster=District, data=Cooper03,
               intervals.type="LB", I2=c("I2q", "I2hm", "I2am", "ICC")) )

## Year as a predictor
summary( meta3(y=y, v=v, cluster=District, x=scale(Year, scale=FALSE), data=Cooper03) )
}
\keyword{datasets}
