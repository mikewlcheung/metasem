\name{Chan17}
\alias{Chan17}
\docType{data}
\title{Dataset from Chan, Jones, Jamieson, and Albarracin (2017)}

\description{A dataset of multiple treatment effects of standardized mean differences on misinformation and debunking effects.}

\usage{data(Chan17)}
\format{
  A data frame with 34 independent samples from 6 research reports.
  \describe{
    \item{\code{Author}}{a character vector of study}
    \item{\code{g_misinfo}}{Hedges' g of misinformation comparing the misinformation experimental and control groups}
    \item{\code{g_debunk}}{Hedges' g of debunking comparing the debuking experimental and misinformation experimental groups}
    \item{\code{v_misinfo}}{sampling variance of g_misinfo}
    \item{\code{c_mis_deb}}{Sampling covariance between \code{g_misinfo} and \code{g_debunk} due to the overlap of the misinformation experimental group}
    \item{\code{v_debunk}}{sampling variance of g_debunk}
    \item{\code{PublicationYear}}{publication year}
    \item{\code{Published}}{published or unpublished}
    \item{\code{MeanAge}}{mean age of participants}
    \item{\code{PctFemale}}{percentage of female participants}
  }
}

\details{
  The sampling variances and covariances are calculated using Gleser and Olkin's (2009) method for multiple treatment effects (Equations 3.3 and 3.4). Since the sample sizes of the misinformation, debunking, and control groups are not given, it is assumed they are equal.}

\source{Chan, M. S., Jones, C. R., Hall Jamieson, K., & Albarracin, D. (2017). Debunking: A meta-analysis of the psychological efficacy of messages countering misinformation. \emph{Psychological Science}, \bold{28(11)}, 1531-1546. https://doi.org/10.1177/0956797617714579}

\references{
Gleser, L. J., & Olkin, I. (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis.} (2nd ed., pp. 357-376). Russell Sage Foundation.
}

\keyword{datasets}
