\name{Hox02}
\alias{Hox02}
\docType{data}
\title{Simulated Effect Sizes Reported by Hox (2002)}
\description{
Twenty stimulated studies on standardized mean difference and one
continuous study characteristic reported by Hox (2002).
}
\usage{data(Hox02)}
\format{
}
\details{
  The variables are:
  \describe{
  \item{study}{Study number}
  \item{yi}{Effect size (standardized mean difference)}
  \item{vi}{Sampling variance of the effect size}
  \item{weeks}{Duration of the experimental intervention in terms of weeks}
}
}
\source{
Hox, J. J. (2002). \emph{Multilevel analysis: Techniques and applications.} Mahwah, N.J.: Lawrence Erlbaum Associates.
}
\references{
Cheung, M. W. L. (2008). A model for integrating fixed-, random-, and mixed-effects meta-analyses into structural equation modeling. \emph{Psychological Methods}, \bold{13}, 182-202. 
}
\examples{
data(Hox02)

#### ML estimation method
## Random-effects meta-analysis
summary( meta(y=yi, v=vi, data=Hox02, I2=c("I2q", "I2hm"), intervals.type="LB") ) 

## Fixed-effects meta-analysis
summary( meta(y=yi, v=vi, data=Hox02, RE.constraints=matrix(0, ncol=1, nrow=1),
              model.name="Fixed effects model") )

## Mixed-effects meta-analysis with "weeks" as a predictor
## Request likelihood-based CI
summary( meta(y=yi, v=vi, x=weeks, data=Hox02, intervals.type="LB",
              model.name="Mixed effects meta analysis with LB CI") )

#### REML estimation method
## Random-effects meta-analysis with REML
summary( VarComp <- reml(y=yi, v=vi, data=Hox02) )

## Extract the variance component
VarComp_REML <- matrix( coef(VarComp), ncol=1, nrow=1 )

## Meta-analysis by treating the variance component as fixed
summary( meta(y=yi, v=vi, data=Hox02, RE.constraints=VarComp_REML) )


## Mixed-effects meta-analysis with "weeks" as a predictor
## Request likelihood-based CI
summary( reml(y=yi, v=vi, x=weeks, intervals.type="LB",
              data=Hox02, model.name="REML with LB CI") )
}
\keyword{datasets}
