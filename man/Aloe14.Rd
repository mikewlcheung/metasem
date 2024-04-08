\name{Aloe14}
\alias{Aloe14}
\docType{data}
\title{Multivariate effect sizes between classroom management
  self-efficacy (CMSE) and other variables reported by Aloe et al. (2014)
}
\description{This study reports sixteen studies on the effect sizes
  (correlation coefficients) between CMSE and emotional exhaustion (EE),
  depersonalization (DP), and (lowered) personal accomplishment (PA)
  reported by Aloe et al. (2014).
}
\usage{data("Aloe14")}
\format{
  A data frame with 16 observations on the following 14 variables.
  \describe{
    \item{\code{Study}}{a factor with levels \code{Betoret} \code{Brouwers & Tomic} \code{Bumen} \code{Chang} \code{Durr} \code{Evers et al.} \code{Friedman} \code{Gold} \code{Huk} \code{Kress} \code{Kumarakulasingam} \code{Martin et al.} \code{Ozdemir} \code{Skaalvik and Skaalvik} \code{Williams}}
    \item{\code{Year}}{Year of publication}
    \item{\code{EE}}{Emotional exhaustion}
    \item{\code{DP}}{Depersonalization}
    \item{\code{PA}}{(Lowered) personal accomplishment}
    \item{\code{V_EE}}{Sampling variance of emotional exhaustion}
    \item{\code{V_DP}}{Sampling variance of depersonalization}
    \item{\code{V_PA}}{Sampling variance of (lowered) personal accomplishment}
    \item{\code{C_EE_DP}}{Sampling covariance between EE and DP}
    \item{\code{C_EE_PA}}{Sampling covariance between EE and PA}
    \item{\code{C_DP_PA}}{Sampling covariance between DP and PA}
    \item{\code{Publication_type}}{Either \code{Dissertation} or \code{Journal}}
    \item{\code{Percentage_females}}{Percentage of females in the study}
    \item{\code{Years_experience}}{Average years of experience}
  }
}
\source{Aloe, A. M., Amo, L. C., & Shanahan, M. E. (2014). Classroom management self-efficacy and burnout: A multivariate meta-analysis. \emph{Educational Psychology Review}, \bold{26(1)}, 101-126. doi:10.1007/s10648-013-9244-0
}
\examples{
\donttest{
data(Aloe14)

## Random-effects meta-analysis
meta1 <- meta(cbind(EE,DP,PA),
              cbind(V_EE, C_EE_DP, C_EE_PA, V_DP, C_DP_PA, V_PA),
              data=Aloe14)
## Remove error code
meta1 <- rerun(meta1)

summary(meta1)

## Extract the coefficients for the variance component of the random effects
coef1 <- coef(meta1, select="random")

## Convert it into a symmetric matrix by row major
my.cov <- vec2symMat(coef1, byrow=TRUE)

## Convert it into a correlation matrix
cov2cor(my.cov)

## Plot the multivariate effect sizes
plot(meta1)
}
}
\keyword{datasets}
