\name{Berkey98}
\alias{Berkey98}
\docType{data}
\title{Five Published Trails from Berkey et al. (1998)}
\description{
  Five published trials, reported by Berkey et al. (1998), comparing surgical and non-surgical treatments
  for medium-severity periodontal disease, one year after treatment. 
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
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
## maybe str(Berkey98) ; plot(Berkey98) ...
}
\keyword{datasets}
