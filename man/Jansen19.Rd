\name{Jansen19}
\alias{Jansen19}
\docType{data}
\title{Effect Sizes Reported by Jansen et al. (2019)}
\description{
  This dataset contains effect sizes derived from Jansen et al. (2019). It includes standardized mean differences (SMDs) regarding self-regulated learning (SRL) activities and achievements, comparing the SRL intervention group to the control group, and the correlations between SRL activities and achievement. There are two key differences compared to the analyses conducted by Jansen et al. (2019). First, while Jansen et al. (2019) used correlation matrices with a d-to-r conversion, this dataset employs both SMDs and correlation matrices. Second, this dataset averages the effect sizes within each study to prevent non-independence among the effect sizes. Consequently, the results from this dataset may not be directly comparable to those of Jansen et al. (2019).
}
\usage{data(Jansen19)}

\details{
  The variables are:
  \describe{
  \item{ArticleID}{Article ID number}
  \item{SRL_act}{SMD of the SRL activity comparing the SRL group to the control group}
  \item{Achieve}{SMD of the achievement comparing the SRL group to the control group}
  \item{r_SRL_achieve}{Correlation between the SRL activity and the achievement}
  \item{Ntotal}{Total sample size of the control and intervention groups}
  \item{Online}{Whether the primary studies were conducted online or offline}
}
}
\source{
Jansen, R. S., van Leeuwen, A., Janssen, J., Jak, S., & Kester, L. (2019). Self-regulated learning partially mediates the effect of self-regulated learning interventions on achievement in higher education: A meta-analysis. \emph{Educational Research Review}, \bold{28}, 100292. https://doi.org/10.1016/j.edurev.2019.100292
}

\examples{
\donttest{
data(Jansen19)
}
}
\keyword{datasets}
