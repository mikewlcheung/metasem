\name{summary}
\alias{summary.tssem1FEM}
\alias{summary.tssem1FEM.cluster}
\alias{summary.tssem1REM}
\alias{summary.wls}
\alias{summary.wls.cluster}
\alias{summary.meta}
\alias{summary.reml}
\alias{print.summary.tssem1FEM}
\alias{print.summary.tssem1FEM.cluster}
\alias{print.summary.wls}
\alias{print.summary.meta}
\alias{print.summary.reml}
\title{Summary Method for tssem1, wls and meta Objects
}
\description{It summaries results for either class \code{tssem1FEM},
  class \code{tssem1FEM.cluster}, class \code{tssem1REM}, class
  \code{wls}, class \code{wls.cluster}, class \code{meta} or class \code{reml}.
}
\usage{
\method{summary}{tssem1FEM}(object, \dots)
\method{summary}{tssem1FEM.cluster}(object, \dots)
\method{summary}{tssem1REM}(object, \dots)
\method{summary}{wls}(object, df.adjustment=0, \dots)
\method{summary}{wls.cluster}(object, df.adjustment=0, \dots)
\method{summary}{meta}(object, homoStat=TRUE, \dots)
\method{summary}{reml}(object, \dots)
\method{print.summary}{tssem1FEM}(x, \dots)
\method{print.summary}{wls}(x, \dots)
\method{print.summary}{meta}(x, \dots)
\method{print.summary}{reml}(x, \dots)
}
\arguments{
  \item{object}{An object returned from either class
	\code{tssem1FEM}, class \code{tssem1FEM.cluster},
	class \code{tssem1REM}, class \code{wls}, class \code{wls.cluster}, class \code{meta} or class \code{reml}.
}
  \item{x}{An object returned from either class \code{summary.tssem1FEM},
	class \code{tssem1FEM.cluster}, class \code{summary.wls}, class \code{summary.meta} or class \code{summary.reml}.
}
  \item{homoStat}{Logical. Whether to conduct a homogeneity test on the
	effect sizes.}
  \item{df.adjustment}{Numeric. Adjust the degrees of freedom
	manually. It may be necessary if the df calculated is incorrect when
	\code{diag.constraints=TRUE}.
	}
  \item{\dots}{Further arguments to be passed to \code{\link[stats]{printCoefmat}}
}
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
%% ~put references to the literature/web site here ~
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
%%  ~~who you are~~
}
\note{If the OpenMx status1 is either 0 or 1, the estimation is
  considered fine. If the OpenMx status1 is other values, it indicates
  estimation problems. Users should refer to \url{http://openmx.psyc.virginia.edu/wiki/errors} for more details.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
  \code{\link[metaSEM]{meta}}, \code{\link[metaSEM]{reml}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

