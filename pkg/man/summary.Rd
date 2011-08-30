\name{summary}
\alias{summary.tssem1FE}
\alias{summary.tssem1FE.cluster}
\alias{summary.tssem1RE}
\alias{summary.wls}
\alias{summary.wls.cluster}
\alias{summary.meta}
\alias{summary.reml}
\alias{print.summary.tssem1FE}
\alias{print.summary.tssem1FE.cluster}
\alias{print.summary.wls}
\alias{print.summary.meta}
\alias{print.summary.reml}
\title{Summary Method for tssem1, wls and meta Objects
}
\description{It summaries results for either class \code{tssem1FE},
  class \code{tssem1FE.cluster}, class \code{tssem1RE}, class
  \code{wls}, class \code{wls.cluster}, class \code{meta} or class \code{reml}.
}
\usage{
\method{summary}{tssem1FE}(object, \dots)
\method{summary}{tssem1FE.cluster}(object, \dots)
\method{summary}{tssem1RE}(object, \dots)
\method{summary}{wls}(object, \dots)
\method{summary}{wls.cluster}(object, \dots)
\method{summary}{meta}(object, \dots)
\method{summary}{reml}(object, \dots)
\method{print.summary}{tssem1FE}(x, \dots)
\method{print.summary}{wls}(x, \dots)
\method{print.summary}{meta}(x, \dots)
\method{print.summary}{reml}(x, \dots)
}
\arguments{
  \item{object}{An object returned from either class
	\code{tssem1FE}, class \code{tssem1FE.cluster},
	class \code{tssem1RE}, class \code{wls}, class \code{wls.cluster}, class \code{meta} or class \code{reml}.
}
  \item{x}{An object returned from either class \code{summary.tssem1FE},
	class \code{tssem1FE.cluster}, class \code{summary.wls}, class \code{summary.meta} or class \code{summary.reml}.
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

\seealso{ \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
  \code{\link[metaSEM]{meta}}, \code{\link[metaSEM]{reml}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

