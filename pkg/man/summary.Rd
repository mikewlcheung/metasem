\name{summary}
\alias{summary.tssem1}
\alias{summary.wls}
\alias{summary.meta}
\alias{print.summary.tssem1}
\alias{print.summary.wls}
\alias{print.summary.meta}
\title{Summary Method for tssem1, wls and meta Objects
}
\description{It summaries results for either \code{tssem1}, \code{wls} or \code{meta} objects.
}
\usage{
\method{summary}{tssem1}(object, \dots)
\method{summary}{wls}(object, \dots)
\method{summary}{meta}(object, \dots)
\method{print.summary}{tssem1}(x, \dots)
\method{print.summary}{wls}(x, \dots)
\method{print.summary}{meta}(x, \dots)
}
\arguments{
  \item{object}{An object returned from either class
	\code{\link[metaSEM]{tssem1}}, class \code{\link[metaSEM]{wls}} or
	class \code{\link[metaSEM]{meta}}
}
  \item{x}{An object returned from either class \code{summary.tssem1}, class \code{summary.wls} or
	class \code{summary.meta}
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}}, \code{\link[metaSEM]{meta}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ summary }

