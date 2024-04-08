\name{Berkey98}
\alias{Berkey98}
\docType{data}
\title{Five Published Trails from Berkey et al. (1998)}
\description{
  The data set includes five published trials, reported by Berkey et al. (1998), comparing surgical and non-surgical treatments
  for medium-severity periodontal disease, one year after treatment. 
}
\usage{data(Berkey98)}

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
\donttest{
data(Berkey98)

#### ML estimation method
## Multivariate meta-analysis
x <- meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL), data=Berkey98)
x <- rerun(x)
summary(x)
plot(x)

## Plot individual studies proportional to the weights
plot(x, study.weight.plot=TRUE)

## Include forest plot from the metafor package
library(metafor)
plot(x, diag.panel=TRUE, main="Multivariate meta-analysis",
     axis.label=c("PD", "AL"))
     forest( rma(yi=PD, vi=var_PD, data=Berkey98) )
     title("Forest plot of PD")
     forest( rma(yi=AL, vi=var_AL, data=Berkey98) )
     title("Forest plot of AL")

## Multivariate meta-analysis with "publication year-1979" as the predictor
summary( meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
              x=scale(pub_year, center=1979), data=Berkey98,
              RE.lbound=NA) )

## Multivariate meta-analysis with equality constraint on the regression coefficients
summary( meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
              x=scale(pub_year, center=1979), data=Berkey98,
              coef.constraints=matrix(c("0.3*Eq_slope", "0.3*Eq_slope"),
              nrow=2)) )

#### REML estimation method
## Multivariate meta-analysis
summary( reml(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
              data=Berkey98,
              model.name="Multivariate meta analysis with REML") )

## Multivariate meta-analysis with "publication year-1979" as the predictor
## Diagonal structure for the variance component
summary( reml(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
              RE.constraints=Diag(c("1e-5*Tau2_1_1", "1e-5*Tau2_2_2")),
              x=scale(pub_year, center=1979), data=Berkey98) )
}
}
\keyword{datasets}
