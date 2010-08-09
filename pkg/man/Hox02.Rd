\name{Hox02}
\alias{Hox02}
\docType{data}
\title{Simulated effect sizes reported by Hox (2002)}
\description{
Twenty stimulated studies on standardized mean difference and one
continuous study characteristic reported by Hox (2002).
}
\usage{data(Hox02)}
\format{
  The format is:
 num [1:20, 1:3] -0.264 -0.23 0.166 0.173 0.225 0.291 0.309 0.435 0.476 0.617 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:3] "yi" "vi" "weeks"
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
## maybe str(Hox02) ; plot(Hox02) ...
}
\keyword{datasets}
