\name{Stadler15}
\alias{Stadler15}
\docType{data}
\title{Correlations from Stadler et al. (2015)}
\description{
  The data set includes correlations between complex problem solving and
  intelligence reported by Stadler et al. (2015). 
}
\usage{data(Stadler15)}

\details{
A list of data with the following structure:
\describe{
  \item{ID}{ID of the effect sizes}
  \item{Authors}{Authors of the studies}
  \item{Year}{Year of the studies}
  \item{N}{Sample size}
  \item{CPSMeasure}{Complex problem solving (CPS) measure}
  \item{IntelligenceMeasure}{Intelligence measure}
  \item{r}{Correlation between CPS and intelligence}
  \item{v}{Sampling variance of r}
}
}
\source{
  Stadler, M., Becker, N., Godker, M., Leutner, D., & Greiff, S. (2015). Complex problem solving and intelligence: A meta-analysis. \emph{Intelligence}, \bold{53}, 92–101.
}

\examples{
data(Stadler15)
}
\keyword{datasets}
