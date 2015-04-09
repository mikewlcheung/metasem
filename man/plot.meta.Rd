\name{plot}
\alias{plot.meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for meta objects
%%  ~~function to do ... ~~
}
\description{It plots the pooled and the individual effect sizes and
  their confidence ellipses for \code{meta} objects.
}
\usage{
\method{plot}{meta}(x, effect.sizes, add.margin = 0.1, interval = 0.95,
main= "Effect Sizes and their Confidence Ellipses",
axis.labels= paste("Effect size ", effect.sizes, sep = ""),
study.col = "black", study.pch = 19, study.min.cex = 0.8,
study.weight.plot = FALSE, study.ellipse.plot = TRUE,
study.ellipse.col = "black", study.ellipse.lty = 2,
study.ellipse.lwd = 0.5,
estimate.col = "blue", estimate.pch = 18, estimate.cex = 2,
estimate.ellipse.plot = TRUE, estimate.ellipse.col = "red",
estimate.ellipse.lty = 1, estimate.ellipse.lwd = 2,
randeff.ellipse.plot = TRUE, randeff.ellipse.col = "green",
randeff.ellipse.lty = 1, randeff.ellipse.lwd = 2,
univariate.plot = TRUE, univariate.lines.col = "gray",
univariate.lines.lty = 3,
univariate.lines.lwd = 1, univariate.polygon.width = 0.02,
univariate.polygon.col = "red", univariate.arrows.col = "green",
univariate.arrows.lwd = 2,
diag.panel = FALSE, xlim=NULL, ylim=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link[metaSEM]{meta}}.
}
  \item{effect.sizes}{Numeric values indicating which effect sizes to
	be plotted. At least two effect sizes are required. To plot the effect sizes of \eqn{y_1}{y1} and
	\eqn{y_2}{y2}, one may use \code{effect.sizes=c(1,2)}. If it is missing, all effect sizes will be plotted in a
	pairwise way.
}
  \item{add.margin}{Value for additional margins on the left and bottom margins.
}
  \item{interval}{Interval for the confidence ellipses.
}
  \item{main}{Main title of each plot. If there are multiple plots, a
	vector of character titles may be used.
}
  \item{axis.labels}{Labels for the effect sizes.
}
\item{study.col}{Color of individual studies. See \code{col} in \code{\link[graphics]{par}}.  
}
  \item{study.pch}{Plotting character of individual studies. See \code{pch} in \code{\link[graphics]{points}}. 
}
  \item{study.min.cex}{Minimum value of cex for individual
	studies. See \code{cex} in \code{\link[graphics]{par}}.
}
  \item{study.weight.plot}{Logical. If \code{TRUE}, the ploting size of individual
	studies (cex) will be proportional to one over the square root of the
	determinant of the sampling covariance matrix of the study.
}
  \item{study.ellipse.plot}{Logical. If \code{TRUE}, the confidence
	ellipses of individual studies are plotted.
}
  \item{study.ellipse.col}{Color of the confidence ellipses of
	individual studies. See \code{col} in \code{\link[graphics]{par}}. 
}
  \item{study.ellipse.lty}{Line type of the confidence ellipse of
	individual studies. See \code{lty} in \code{\link[graphics]{par}}. 
}
  \item{study.ellipse.lwd}{Line width of the confidence ellipse of
	individual studies. See \code{lwd} in \code{\link[graphics]{par}}. 
}
  \item{estimate.col}{Color of the estimated effect size. See \code{col} in \code{\link[graphics]{par}}. 
}
  \item{estimate.pch}{Plotting character of the estimated effect sizes. See \code{pch} in \code{\link[graphics]{points}}.
}
  \item{estimate.cex}{Amount of plotting of the estimated effect sizes. See \code{cex} in \code{\link[graphics]{par}}.
}
  \item{estimate.ellipse.plot}{Logical. If \code{TRUE}, the confidence
	ellipse of the estimated effect sizes will be plotted.
}
  \item{estimate.ellipse.col}{Color of the confidence
	ellipse of the estimated effect sizes. See \code{col} in \code{\link[graphics]{par}}. 
}
  \item{estimate.ellipse.lty}{Line type of the confidence
	ellipse of the estimated effect sizes. See \code{lty} in \code{\link[graphics]{par}}. 
}
  \item{estimate.ellipse.lwd}{Line width of the confidence
	ellipse of the estimated effect sizes. See \code{lwd} in \code{\link[graphics]{par}}. 
}
  \item{randeff.ellipse.plot}{Logical. If \code{TRUE}, the confidence
	ellipses of the random effects will be plotted.
}
  \item{randeff.ellipse.col}{Color of the confidence
	ellipses of the random effects. See \code{col} in \code{\link[graphics]{par}}. 
}
  \item{randeff.ellipse.lty}{Line type of the confidence
	ellipses of the random effects. See \code{lty} in \code{\link[graphics]{par}}. 
}
  \item{randeff.ellipse.lwd}{Line width of the confidence
	ellipses of the random effects. See \code{lwd} in \code{\link[graphics]{par}}. 
}
  \item{univariate.plot}{Logical. If \code{TRUE}, the estimated univariate effect
	sizes will be plotted.
}
  \item{univariate.lines.col}{Color of the estimated univariate effect
	sizes. See \code{col} in \code{\link[graphics]{par}}. 
}
  \item{univariate.lines.lty}{Line type of the estimated univariate effect
	sizes. See \code{lty} in \code{\link[graphics]{par}}. 
}
  \item{univariate.lines.lwd}{Line width of the estimated univariate effect
	sizes. See \code{lwd} in \code{\link[graphics]{par}}. 
}
  \item{univariate.polygon.width}{Width of the polygon of the estimated univariate effect
	sizes.
}
  \item{univariate.polygon.col}{Color of the polygon of the estimated univariate effect
	sizes.
}
  \item{univariate.arrows.col}{Color of the arrows of the estimated univariate effect
	sizes.
}
  \item{univariate.arrows.lwd}{Line width of the arrows of the estimated univariate effect
	sizes.
}
  \item{diag.panel}{Logical. If \code{TRUE}, diagonal panels will be
	created. They can then be used for forrest plots for univariate meta-analysis.
}
  \item{xlim}{NULL or a numeric vector of length 2; if it is NULL, it
	provides defaults estimated from the data.}
  \item{ylim}{NULL or a numeric vector of length 2; if it is NULL, it
	provides defaults estimated from the data.}
  \item{\dots}{Further arguments to be passed to the plot function.
}
}


\references{
  Cheung, M. W.-L. (2013). Multivariate meta-analysis as structural equation models. \emph{Structural Equation Modeling}, \bold{20}, 429-454.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{The estimated effect sizes and random effects are based on the
  labels Intercept1, Intercept2, ... and Tau2_1_1, Tau2_2_1, Tau2_2_2,
  etc. At least two effect sizes are required for this function. The
  .ellipse() function is adopted from the ellipse package written by
  Duncan Murdoch and E. D. Chow.  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\dontrun{
## Multivariate meta-analysis
x <- meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL), data=Berkey98)
plot(x)

## Plot individual studies proportional to the weights
plot(x, study.weight.plot=TRUE)

## Include forest plot from the metafor package
library(metafor)
plot(x, diag.panel=TRUE, main="Multivariate meta-analysis",
axis.label=c("PD", "AL"))
forest( rma(yi=PD, vi=var_PD, data=Berkey98) )
title("Forest plot of PD")
forest( rma(yi=AL, vi=var_AL, data=Berkey98) )
title("Forest plot of AL")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }% __ONLY ONE__ keyword per line
