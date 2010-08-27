\name{Digman97}
\alias{Digman97}
\alias{Digman97.n}
\alias{Digman97.cluster}
\docType{data}
\title{Factor correlation matrices of the Big Five Model from Digman (1997)}
\description{
Fourteen studies of the factor correlation matrices of the Big Five Model of personality reported by Digman (1997). 
}
\usage{data(Digman97)}
\format{
  The format is:
List of 14
 $ Digman 1 (1994)                : num [1:5, 1:5] 1 -0.48 -0.1 0.27 0.37 -0.48 1 0.62 0.41 0 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:5] "E" "A" "C" "ES" ...
  .. ..$ : chr [1:5] "E" "A" "C" "ES" ...
 $ Digman 2 (1994)                : num [1:5, 1:5] 1 -0.3 0.07 0.09 0.45 -0.3 1 0.39 0.53 -0.05 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:5] "E" "A" "C" "ES" ...
  .. ..$ : chr [1:5] "E" "A" "C" "ES" ...
...
}
\details{
The five factors are:
\describe{
  \item{E}{Extraversion}
  \item{A}{Agreeableness}
  \item{C}{Conscientiousness}
  \item{ES}{Emotional Stability}
  \item{I}{Intellect or Openness to Experience}
}
  Digman97.n. and Digman97.cluster store the sample sizes and the types
  of participants of the studies, respectively.
}
\source{
Digman, J.M. (1997). Higher-order factors of the Big Five. \emph{Journal of Personality and Social Psychology}, \bold{73}, 1246-1256.
}
\references{
 Cheung, M. W. L., & Chan, W. (2005). Classifying correlation matrices into relatively homogeneous subgroups: A cluster analytic approach. \emph{Educational and Psychological Measurement}, \bold{65}, 954-979.
}
\examples{
data(Digman97)
## maybe str(Digman97)  ...
}
\keyword{datasets}
