\name{issp05}
\alias{issp05}
\docType{data}
\title{A Dataset from ISSP (2005) 
}
\description{Thirty-two covariance matrices on work-related attitudes
  were extracted from the International Social Survey Programme 2005:
  Work Orientation III (ISSP, 2005). Seven variables were selected for
  demonstration purposes. They were grouped into three constructs:
  \emph{Importance of Job Prospects} measured by job security (JP1),
  high income (JP2), and opportunity for advancement (JP3); \emph{Importance
  of Job Autonomy} measured by work independently (JA1) and decide time
  of work (JA2); and \emph{Importance of Contributions to Society} measured by help other
  people (CS1) and a job useful to society (CS2).
}
\usage{data(issp05)}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 32 covariance matrices }
  \item{n}{A vector of sample sizes}
  \item{means}{A matrix of means}
  \item{pdi}{Hofstede's Power Distance Index}
  \item{idv}{Hofstede's Individualism}
  \item{mas}{Hofstede's Masculinity}
  \item{uai}{Hofstede's Uncertainty Avoidance Index}
  \item{ltowvs}{Hofstede's Long- Versus Short-Term Orientation}
  \item{ivr}{Hofstede's Indulgence Versus Restraint}
}
}
\source{ISSP Research Group (2007): \emph{International Social Survey Programme 2005: Work Orientation III (ISSP 2005)}. GESIS Data Archive, Cologne. ZA4350 Data file Version 1.0.0, doi:10.4232/1.4350

 https://geerthofstede.com/research-and-vsm/dimension-data-matrix/
}


\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{issp89}}
}

\examples{

\donttest{
data(issp05)

#### TSSEM random-effects model with covariance matrices

## Stage 1 analysis
rand1 <- tssem1(issp05$data, issp05$n, method="REM", cor.analysis=FALSE)
summary(rand1)

## Proposed model
model1 <- "JP =~ JP1 + JP2 + JP3
           JA =~ JA1 + JA2
           CS =~ CS1 + CS2
           JP ~~ JA + CS
           JA ~~ CS"

ram1 <- lavaan2RAM(model1, obs.variables=c("JP1", "JP2", "JP3", "JA1", "JA2",
                                           "CS1", "CS2"))

## Stage 2 analysis
rand2 <- tssem2(rand1, RAM=ram1)
summary(rand2)

plot(rand2)

#### OSMASEM with covariance matrices
## Create a data frame for the OSMASEM
df <- Cor2DataFrame(issp05$data, n=issp05$n, Means=issp05$means,
                    cor.analysis=FALSE)

## Standardize idv
idv <- scale(issp05$idv)

## Replace missing values with mean
idv[is.na(idv)] <- mean(idv, na.rm=TRUE)
df$data$idv <- idv

## No moderator 
fit1 <- osmasem2(model.name="No_moderator", RAM=ram1, data=df,
                 cor.analysis=FALSE, mean.analysis=FALSE)
summary(fit1, fitIndices = TRUE)

## Proposed model with idv as a moderator
model2 <- "JP =~ a*JP1 + b*JP2 + c*JP3
           JA =~ d*JA1 + e*JA2
           CS =~ f*CS1 + g*CS2
           JP ~~ JA + CS
           JA ~~ CS
           a == a0 + a1*data.idv
           b == b0 + b1*data.idv
           c == c0 + c1*data.idv
           d == d0 + d1*data.idv
           e == e0 + e1*data.idv
           f == f0 + f1*data.idv
           g == g0 + g1*data.idv"

ram2 <- lavaan2RAM(model2, obs.variables=c("JP1", "JP2", "JP3", "JA1", "JA2",
                                           "CS1", "CS2"))

fit2 <- osmasem2(RAM=ram2, data=df, cor.analysis=FALSE, mean.analysis=FALSE,
                 replace.constraints = TRUE)
summary(fit2)

## Compare fit1 and fit2
anova(fit2, fit1)
}
}
\keyword{datasets}
