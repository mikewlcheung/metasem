\name{rerun}
\alias{rerun}
\title{Rerun models via mxRun()
}
\description{It reruns models via mxRun().
}
\usage{
rerun(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of either class \code{tssem1FEM},
  class \code{tssem1REM}, class \code{wls}, class \code{meta} or class \code{reml}.
}
  \item{\dots}{{Futher arguments to be passed to \code{\link[OpenMx]{mxRun}}}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
random1 <- tssem1(Digman97$data, Digman97$n, method="REM", RE.type="Diag")
random1_rerun <- rerun(random1)
summary(random1_rerun)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{meta-analysis}
\keyword{tssem}
