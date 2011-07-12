\name{wvs94}
\alias{wvs94}
\docType{data}
\title{Forty-four Studies from Cheung (2011)}
\description{
Between 1990 and 1993, 57,561 adults aged 18 and above from 42 nations
were interviewed by local academic institutes in Eastern European nations and by professional
survey organizations in other nations.
Standardized mean difference (SMD) between males and females on life
satisfaction and on life control in each country were calculated as the effect sizes. Positive values indicate that males have higher scores than females do.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(wvs94)}
\format{
}
\details{
    The variables are:
  \describe{
  \item{country}{Country}
  \item{lifesat}{SMD on life satisfaction}
  \item{lifecon}{SMD on life control}
  \item{lifesat.var}{Sampling variance of lifesat}
  \item{inter.cov}{Sampling covariance between lifesat and lifecon}
  \item{lifecon.var}{Sampling variance of lifecon}
  \item{gnp}{Gross National Product}
}

}
\source{
  World Values Study Group. (1994). World Values Survey, 1981-1984 and 1990-1993 [Computer file]. \emph{Ann Arbor, MI: Inter-university Consortium for Political and Social Research.}
}
\references{
  Au, K., & Cheung, M. W. L. (2004). Intra-cultural variation and job autonomy in 42
  countries. \emph{Organization Studies}, \bold{25}, 1339-1362.
  
  Cheung, M. W. L. (2011). Multivariate meta-analaysis as structural
  equation models. \emph{Manuscript submitted for publication.}
}
\examples{
data(wvs94)
## maybe str(wvs94) ; plot(wvs94) ...
}
\keyword{datasets}
