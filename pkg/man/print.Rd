\name{print}
\alias{print.tssem1FE}
\alias{print.tssem1FE.cluster}
\alias{print.tssem1RE}
\alias{print.wls}
\alias{print.meta}
\alias{print.reml}
\title{Print Methods for \code{tssem1FE}, \code{tssem1FE.cluster}, \code{tssem1RE}, \code{wls}, \code{meta} and \code{reml} Objects 
}
\description{Print methods for either \code{tssem1FE},
  \code{tssem1FE.cluster}, \code{tssem1RE}, \code{wls},
  \code{meta} or \code{reml} objects.
}
\usage{
\method{print}{tssem1FE}(x, \dots)
\method{print}{tssem1FE.cluster}(x, \dots)
\method{print}{tssem1RE}(x, \dots)
\method{print}{wls}(x, \dots)
\method{print}{meta}(x, \dots)
\method{print}{reml}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object returned from either class
	\code{tssem1FE}, class \code{tssem1FE.cluster}, class \code{tssem1RE},
	class \code{wls}, class \code{meta} or class \code{reml}
}
  \item{\dots}{Further arguments to be passed to \code{summary.default}
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
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
  \code{\link[metaSEM]{meta}},  \code{\link[metaSEM]{reml}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ methods }
