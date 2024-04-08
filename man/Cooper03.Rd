\name{Cooper03}
\alias{Cooper03}
\docType{data}
\title{Selected effect sizes from Cooper et al. (2003)}
\description{
  Fifty-six effect sizes from 11 districts from Cooper et al. (2003) were reported by Konstantopoulos (2011).
}
\usage{data(Cooper03)}
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
\donttest{
data(Cooper03)

#### ML estimation method
## No predictor
summary( model1 <- meta3L(y=y, v=v, cluster=District, data=Cooper03) )

## Show all heterogeneity indices and their 95\% confidence intervals
summary( meta3L(y=y, v=v, cluster=District, data=Cooper03,
               intervals.type="LB", I2=c("I2q", "I2hm", "I2am", "ICC")) )

## Year as a predictor
summary( meta3L(y=y, v=v, cluster=District, x=scale(Year, scale=FALSE),
               data=Cooper03, model.name="Year as a predictor") )

## Equality of level-2 and level-3 heterogeneity
summary( model2 <- meta3L(y=y, v=v, cluster=District, data=Cooper03,
                         RE2.constraints="0.2*EqTau2",
                         RE3.constraints="0.2*EqTau2",
                         model.name="Equal Tau2") )

## Compare model2 vs. model1
anova(model1, model2)

#### REML estimation method
## No predictor
summary( reml3L(y=y, v=v, cluster=District, data=Cooper03) )

## Level-2 and level-3 variances are constrained equally 
summary( reml3L(y=y, v=v, cluster=District, data=Cooper03,
               RE.equal=TRUE, model.name="Equal Tau2") )

## Year as a predictor
summary( reml3L(y=y, v=v, cluster=District, x=scale(Year, scale=FALSE),
               data=Cooper03, intervals.type="LB") )

## Handling missing covariates with FIML
## Create 20/56 MCAR data in Year
set.seed(10000)
Year_MCAR <- Cooper03$Year
Year_MCAR[sample(56, 20)] <- NA
summary( meta3LFIML(y=y, v=v, cluster=District, x2=scale(Year_MCAR, scale=FALSE),
                    data=Cooper03, model.name="NA in Year_MCAR") )
}
}
\keyword{datasets}
