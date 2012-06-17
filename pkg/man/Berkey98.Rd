\name{Berkey98}
\alias{Berkey98}
\docType{data}
\title{Five Published Trails from Berkey et al. (1998)}
\description{
  The data set includes five published trials, reported by Berkey et al. (1998), comparing surgical and non-surgical treatments
  for medium-severity periodontal disease, one year after treatment. 
}
\usage{data(Berkey98)}
\format{
}
\details{
    The variables are:
  \describe{
  \item{trial}{Trial number}
  \item{pub_year}{Publication year}
  \item{no_of_patients}{Number of patients}
  \item{PD}{Patient improvements (mm) in \emph{probing depth}}
  \item{AL}{Patient improvements (mm) in \emph{attachment level}}
  \item{var_PD}{Sampling variance of PD}
  \item{cov_PD_AL}{Sampling covariance between PD and AD}
  \item{var_AL}{Sampling variance of AL}
}

}
\source{
  Berkey, C. S., Hoaglin, D. C., Antczak-Bouckoms, A., Mosteller, F, & Colditz, G. A. (1998). Meta-analysis of multiple outcomes by regression with random effects. \emph{Statistics in Medicine}, \bold{17}, 2537-2550.
}
\examples{
data(Berkey98)

#### ML estimation method
## Multivariate meta-analysis
summary( meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
              data=Berkey98, I2=c("I2q", "I2hm", "I2am"),
              intervals.type="LB") )

## Multivariate meta-analysis with "publication year-1979" as the predictor
summary( meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
              x=scale(pub_year, center=1979), data=Berkey98) )

## Multivariate meta-analysis with equality constraint on the regression coefficients
summary( meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
              x=scale(pub_year, center=1979), data=Berkey98,
              coef.constraints=matrix(c("0.3*Eq_slope", "0.3*Eq_slope"),
              nrow=2)) )

#### REML estimation method
## Multivariate meta-analysis
summary( reml(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL,var_AL),
              data=Berkey98,
              model.name="Multivariate meta analysis with REML") )

## Multivariate meta-analysis with "publication year-1979" as the predictor
summary( reml(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
              x=scale(pub_year, center=1979), data=Berkey98) )
}
\keyword{datasets}
