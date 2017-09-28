\name{Nohe15}
\alias{Nohe15A1}
\alias{Nohe15A2}
\docType{data}
\title{Correlation Matrices from Nohe et al. (2015)}
\description{
The data sets include two lists of correlation matrices of panel studies between
work-family conflict and strain reported in Table A1 (\code{Nohe15A1})
and Table A2 (\code{Nohe15A2}) by Nohe et al. (2015). 
}
\usage{data(Nohe15A1)
data(Nohe15A2)}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of studies of correlation matrices. The
	variables are \emph{W1}, \emph{S1}, \emph{W2}, and \emph{S2} in
	\code{Nohe15A1} and \emph{F1}, \emph{S1}, \emph{F2}, and \emph{S2} in
	\code{Nohe15A2}}
  \item{n}{A vector of sample sizes}
  \item{RelXX}{The reliabilities of \emph{W1}, \emph{S1}, \emph{W2} and
	\emph{S2} in \code{Nohe15A1} and the reliabilities of \emph{F1}
	\emph{S1}, \emph{F2} , and \emph{S2} in \code{Nohe15A2}}
  \item{FemalePer}{Percentage of female participants}
  \item{Publication}{Whether the studies were published (\emph{P}) or
	unpublished (\emph{U})}
  \item{Lag}{Time lag between the coded measurement waves in months}
}
}
\source{
  Nohe, C., Meier, L. L., Sonntag, K., & Michel, A. (2015). The chicken or the egg? A meta-analysis of panel studies of the relationship between work-family conflict and strain. \emph{Journal of Applied Psychology}, \bold{100}(2), 522-536.
}
\examples{
\dontrun{
## Table A1
randA1a <- tssem1(Nohe15A1$data, Nohe15A1$n, method="REM", RE.type="Diag")
summary(randA1a)

model1 <- 'W2 ~ w2w*W1 + s2w*S1
           S2 ~ w2s*W1 + s2s*S1
           W1 ~~ w1WITHs1*S1
           W2 ~~ w2WITHs2*S2
           W1 ~~ 1*W1
           S1 ~~ 1*S1
           W2 ~~ Errw2*W2
           S2 ~~ Errs2*S2'

RAM1 <- lavaan2RAM(model1, obs.variables=c("W1", "S1", "W2", "S2"))
RAM1

randA1b <- tssem2(randA1a, Amatrix=RAM1$A, Smatrix=RAM1$S)
summary(randA1b)

my.plotA1 <- meta2semPlot(randA1b)
semPlot::semPaths(my.plotA1, whatLabels="est", sizeMan=8, edge.label.cex=1.5, color="yellow",
                  edge.color = "black", weighted=FALSE)

## Table A2
randA2a <- tssem1(Nohe15A2$data, Nohe15A2$n, method="REM", RE.type="Diag")
## Rerun to remove error code
randA2a <- rerun(randA2a)
summary(randA2a)

model2 <- 'F2 ~ f2f*F1 + s2F*S1
           S2 ~ f2s*F1 + s2s*S1
           F1 ~~ f1WITHs1*S1
           F2 ~~ f2WITHs2*S2
           F1 ~~ 1*F1
           S1 ~~ 1*S1
           F2 ~~ Errf2*F2
           S2 ~~ Errs2*S2'

RAM2 <- lavaan2RAM(model2, obs.variables=c("F1", "S1", "F2", "S2"))
RAM2

randA2b <- tssem2(randA2a, Amatrix=RAM2$A, Smatrix=RAM2$S)
summary(randA2b)

my.plotA2 <- meta2semPlot(randA2b)
semPlot::semPaths(my.plotA2, whatLabels="est", sizeMan=8, edge.label.cex=1.5, color="yellow",
                  edge.color = "black", weighted=FALSE)

## Parametric bootstrap based on Yu et al. (2016)
## I assume that you know what you are doing!

## Set seed for reproducibility
set.seed(39128482)

## Average the correlation coefficients with the univariate-r approach
uni1 <- uniR1(Nohe15A1$data, Nohe15A1$n)
uni1

## Generate random correlation matrices
boot.cor <- bootMASEM1(uni1, Rep=50)

## Display the quality of the generated correlation matrices
summary(boot.cor)

## Proposed saturated model
model1 <- 'W2 + S2 ~ W1 + S1'

## Use the harmonic mean of the sample sizes as n in SEM
n <- uni1$n.harmonic    
    
boot.fit1 <- bootMASEM2(model=model1, data=boot.cor, n=n)
summary(boot.fit1)

## Proposed model with equal regression coefficients
model2 <- 'W2 ~ Same*W1 + Cross*S1
           S2 ~ Cross*W1 + Same*S1'

boot.fit2 <- bootMASEM2(model=model2, data=boot.cor, n=n)
summary(boot.fit2)
    
}
}
\keyword{datasets}
