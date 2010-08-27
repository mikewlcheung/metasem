\name{Cheung09}
\alias{Cheung09}
\alias{Cheung09.n}
\docType{data}
\title{Data set from TSSEM User's Guide version 1.11 by Cheung (2009) 
}
\description{Four studies were selected from the data set used by Cheung
  and Chan (2005; 2009). Some variables were randomly deleted to
  illustrate the analysis with missing data.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(Cheung09)}
\format{
  The format is:
List of 4
 $ 1: num [1:9, 1:9] 0.773 0.27 0.24 0.238 0.209 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:9] "x1" "x2" "x3" "x4" ...
  .. ..$ : chr [1:9] "x1" "x2" "x3" "x4" ...
 $ 2: num [1:9, 1:9] NA NA NA NA NA NA NA NA NA NA ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:9] "x1" "x2" "x3" "x4" ...
  .. ..$ : chr [1:9] "x1" "x2" "x3" "x4" ...
 $ 3: num [1:9, 1:9] 0.958 0.33 NA NA NA ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:9] "x1" "x2" "x3" "x4" ...
  .. ..$ : chr [1:9] "x1" "x2" "x3" "x4" ...
 $ 4: num [1:9, 1:9] 0.84 0.211 0.142 0.133 0.179 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:9] "x1" "x2" "x3" "x4" ...
  .. ..$ : chr [1:9] "x1" "x2" "x3" "x4" ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{Cheung, M. W. L. (2009). TSSEM: A LISREL syntax generator for two-stage structural equation modeling (Version 1.11) [Computer software]. Retrieved from \url{http://courses.nus.edu.sg/course/psycwlm/internet/tssem.zip}.
}
\references{
  Cheung, M. W. L., & Chan, W. (2005). Meta-analytic structural equation modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10}, 40-64.

Cheung, M. W. L., & Chan, W. (2009). A two-stage approach to synthesizing covariance matrices in meta-analytic structural equation modeling. \emph{Structural Equation Modeling}, \bold{16}, 28-53.
}
\examples{
data(Cheung09)
}
\keyword{datasets}
