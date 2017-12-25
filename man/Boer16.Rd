\name{Boer16}
\alias{Boer16}
\docType{data}
\title{Correlation Matrices from Boer et al. (2016)}
\description{
  The data set includes correlation matrices of leader-member exchange in
  transformational leadership reported by Boer et al. (2016). 
}
\usage{data(Boer16)}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of correlation matrices. The variables are
	\emph{LMX} (leader-member exchange), \emph{TFL}
	(transformational leadership), \emph{JS} (job satisfaction),
	\emph{OC} (organizational commitment), and \emph{LE} (leader
	effectiveness)}
  \item{n}{A vector of sample sizes}
  \item{RelLMX}{The reliability of \emph{LMX}}
  \item{RelTFL}{The reliability of \emph{TFL}}
}
}
\source{
  Boer, D., Deinert, A., Homan, A. C., & Voelpel, S. C. (2016). Revisiting the mediating role of leader-member exchange in transformational leadership: the differential impact model. \emph{European Journal of Work and Organizational Psychology}, \bold{25}(6), 883-899.
}
\examples{
\dontrun{
## Stage 1 analysis
rand1 <- tssem1(Boer16$data, Boer16$n, method="REM", RE.type="Diag", acov="individual")
summary(rand1)

## Stage 2 analysis
model2a <- 'JS+OC+LE ~ LMX+TFL
            LMX ~ TFL
            ## Variance of TFL is fixed at 1
            TFL ~~ 1*TFL
            ## Correlated residuals
            JS ~~ OC
            JS ~~ LE
            OC ~~ LE'

RAM2a <- lavaan2RAM(model2a, obs.variables = c("LMX", "TFL", "JS", "OC", "LE"),
                    A.notation="on", S.notation="with")

rand2a <- tssem2(rand1, Amatrix=RAM2a$A, Smatrix=RAM2a$S)
summary(rand2a)

my.plot2a <- meta2semPlot(rand2a)
semPlot::semPaths(my.plot2a, whatLabels="est", sizeMan=8,
                  edge.label.cex=1.5, color="yellow",
                  edge.color = "black", weighted=FALSE, layout="tree2")
}
}
\keyword{datasets}
