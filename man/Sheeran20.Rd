\name{Sheeran20}
\alias{Sheeran20}
\docType{data}
\title{Effect Sizes Reported by Sheeran20 et al. (2020)}
\description{
  This dataset contains effect sizes derived from Sheeran et al. (2020). It includes standardized mean differences (SMDs) regarding health behaviors, autonomous motivation, and perceived competence, comparing the self-determination theory (SDT) group and control group, and the correlations among the variables. One key difference from the analyses conducted by Sheeran et al. (2020) is that while they used correlation matrices with a d-to-r conversion, this dataset utilizes both standardized mean differences (SMDs) and correlation matrices.
}
\usage{data(Sheeran20)}

\details{
  The variables are:
  \describe{
    \item{Study}{Study ID number}
    \item{n_c}{Cluster-adjusted N in the control group}
    \item{n_t}{Cluster-adjusted N in the treatment group}
    \item{Heal_beh}{SMD of the health behaviors between the treatmeant and control groups}
    \item{Auto_mot}{SMD of the autonomous motivation between the treatment and control groups}
    \item{Per_com}{SMD of the perceived competence between the treatment and control groups}
    \item{r_Auto_mot_Heal_beh}{correlation between the autonomous motivation and the health behaviors}
    \item{r_Per_com_Heal_beh}{correlation between the perceived competence and the health behaviors}
    \item{r_Per_com_Auto_mot}{correlation between the perceived competence and the autonomous motivation}
    \item{Beh_timing}{Time between intervention starting and measurement of behavior in weeks}
    \item{SDT_timing}{Time between intervention starting and post-intervention measurement of SDT constructs in weeks}
    \item{Overweight}{Overweight/obese participants composed entire sample}
    \item{College}{University/college students composed entire sample}
    \item{Adolescents}{Adolescents composed entire sample}
    \item{Clinical_sample}{Older adults composed entire sample}
}
}
\source{
  Sheeran, P., Wright, C. E., Avishai, A., Villegas, M. E., Lindemans, J. W., Klein, W. M. P., Rothman, A. J., Miles, E., & Ntoumanis, N. (2020). Self-determination theory interventions for health behavior change: Meta-analysis and meta-analytic structural equation modeling of randomized controlled trials. \emph{Journal of Consulting and Clinical Psychology}, \bold{88}(8), 726-737. https://doi.org/10.1037/ccp0000501
}

\examples{
\donttest{
data(Sheeran20)
}
}
\keyword{datasets}
