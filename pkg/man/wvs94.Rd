\name{wvs94}
\alias{wvs94}
\docType{data}
\title{Forty-four Studies from Cheung (in press)}
\description{
Between 1990 and 1993, 57,561 adults aged 18 and above from 42 nations
were interviewed by local academic institutes in Eastern European nations and by professional
survey organizations in other nations.
Standardized mean difference (SMD) between males and females on life
satisfaction and on life control in each country were calculated as the effect sizes. Positive values indicate that males have higher scores than females do.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(wvs94)}
\format{
}
\details{
    The variables are:
  \describe{
  \item{country}{Country}
  \item{lifesat}{SMD on life satisfaction}
  \item{lifecon}{SMD on life control}
  \item{lifesat.var}{Sampling variance of lifesat}
  \item{inter.cov}{Sampling covariance between lifesat and lifecon}
  \item{lifecon.var}{Sampling variance of lifecon}
  \item{gnp}{Gross National Product}
}

}
\source{
  World Values Study Group. (1994). World Values Survey, 1981-1984 and 1990-1993 [Computer file]. \emph{Ann Arbor, MI: Inter-university Consortium for Political and Social Research.}
}
\references{
  Au, K., & Cheung, M. W. L. (2004). Intra-cultural variation and job autonomy in 42
  countries. \emph{Organization Studies}, \bold{25}, 1339-1362.
  
  Cheung, M. W. L. (in press). Multivariate meta-analaysis as structural
  equation models. \emph{Structural Equation Modeling.}
}
\examples{
data(wvs94)

## Random-effects model
random.ma1 <- meta(y=cbind(lifesat, lifecon),
                   v=cbind(lifesat.var, inter.cov, lifecon.var), data=wvs94,
                   model.name="Random effects model")
summary(random.ma1)

## Random-effects model with both population effect sizes at 0
random.ma2 <- meta(y=cbind(lifesat, lifecon),
                   v=cbind(lifesat.var, inter.cov, lifecon.var), data=wvs94,
                   intercept.constraints=matrix(0, nrow=1, ncol=2),
                   model.name="Effect sizes are fixed at 0")
summary(random.ma2)

## Compare the nested models
anova(random.ma1, random.ma2)

## Fixed-effects model by fixing the variance component at 0 
fixed.ma <- meta(y=cbind(lifesat, lifecon),
                 v=cbind(lifesat.var, inter.cov, lifecon.var), data=wvs94,
                 RE.constraints=matrix(0, ncol=2, nrow=2),
                 model.name="Fixed effects model")
summary(fixed.ma)


## Mixed-effects model
## gnp is divided by 10000 and centered by using 
## scale(gnp/10000, scale=FALSE)
mixed.ma1 <- meta(y=cbind(lifesat, lifecon),
                  v=cbind(lifesat.var, inter.cov, lifecon.var),
                  x=scale(gnp/10000, scale=FALSE), data=wvs94,
                  model.name="GNP as a predictor")
summary(mixed.ma1)

## Mixed-effects model with equal regression coefficients
mixed.ma2 <- meta(y=cbind(lifesat, lifecon),
                  v=cbind(lifesat.var, inter.cov, lifecon.var),
                  x=scale(gnp/10000, scale=FALSE), data=wvs94,
                  coeff.constraints=matrix(c("0.0*Eq_slope",
                                             "0.0*Eq_slope"), nrow=2),
                  model.name="GNP as a predictor with equal slope")
summary(mixed.ma2)

## Compare the nested models
anova(mixed.ma1, mixed.ma2)

## Plot the multivariate effect sizes
plot(random.ma1, main="Estimated effect sizes and their 95\% confidence ellipses",
     axis.label=c("Gender difference on life satisfaction",
                  "Gender difference on life control"))

}
\keyword{datasets}
