\name{wvs94b}
\alias{wvs94b}
\docType{data}
\title{Forty-four Covariance Matrices on Life Satisfaction, Job Satisfaction, and Job Autonomy}
\description{
Between 1990 and 1993, 57,561 adults aged 18 and above from 42 nations
were interviewed by local academic institutes in Eastern European
nations and by professional survey organizations in other nations. The
covariance matrices among Life Satisfaction, Job Satisfaction, and Job Autonomy were calculated.
}
\usage{data(wvs94b)}

\details{
    The variables are:
  \describe{
	\item{data}{Covariance matrix among Life Satisfaction (LS), Job
    Satisfaction (JS), and Job Autonomy (JA)}
  \item{n}{Sample size in the country}
}

}
\source{
  World Values Study Group. (1994). World Values Survey, 1981-1984 and 1990-1993 [Computer file]. \emph{Ann Arbor, MI: Inter-university Consortium for Political and Social Research.}
}
\references{
  Au, K., & Cheung, M. W.-L. (2004). Intra-cultural variation and job autonomy in 42
  countries. \emph{Organization Studies}, \bold{25}, 1339-1362.
}
\examples{
data(wvs94b)

## Get the indirect and the direct effects and
## their sampling covariance matrices for each study
indirect1 <- indirectEffect(wvs94b$data, wvs94b$n)
indirect1

## Multivariate meta-analysis on the indirect and direct effects
indirect2 <- meta(indirect1[, c("ind_eff", "dir_eff")],
                  indirect1[, c("ind_var", "ind_dir_cov", "dir_var")])

summary(indirect2)  
}
\keyword{datasets}
