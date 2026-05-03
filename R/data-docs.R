

#' Multivariate effect sizes between classroom management self-efficacy (CMSE)
#' and other variables reported by Aloe et al. (2014)
#'
#' This study reports sixteen studies on the effect sizes (correlation
#' coefficients) between CMSE and emotional exhaustion (EE), depersonalization
#' (DP), and (lowered) personal accomplishment (PA) reported by Aloe et al.
#' (2014).
#'
#'
#' @name Aloe14
#' @docType data
#' @format A data frame with 16 observations on the following 14 variables.
#' \describe{ \item{\code{Study}}{a factor with levels \code{Betoret}
#' \code{Brouwers & Tomic} \code{Bumen} \code{Chang} \code{Durr} \code{Evers et
#' al.} \code{Friedman} \code{Gold} \code{Huk} \code{Kress}
#' \code{Kumarakulasingam} \code{Martin et al.} \code{Ozdemir} \code{Skaalvik
#' and Skaalvik} \code{Williams}} \item{\code{Year}}{Year of publication}
#' \item{\code{EE}}{Emotional exhaustion} \item{\code{DP}}{Depersonalization}
#' \item{\code{PA}}{(Lowered) personal accomplishment}
#' \item{\code{V_EE}}{Sampling variance of emotional exhaustion}
#' \item{\code{V_DP}}{Sampling variance of depersonalization}
#' \item{\code{V_PA}}{Sampling variance of (lowered) personal accomplishment}
#' \item{\code{C_EE_DP}}{Sampling covariance between EE and DP}
#' \item{\code{C_EE_PA}}{Sampling covariance between EE and PA}
#' \item{\code{C_DP_PA}}{Sampling covariance between DP and PA}
#' \item{\code{Publication_type}}{Either \code{Dissertation} or
#' \code{Journal}} \item{\code{Percentage_females}}{Percentage of females in
#' the study} \item{\code{Years_experience}}{Average years of experience} }
#' @source Aloe, A. M., Amo, L. C., & Shanahan, M. E. (2014). Classroom
#' management self-efficacy and burnout: A multivariate meta-analysis.
#' \emph{Educational Psychology Review}, \bold{26(1)}, 101-126.
#' doi:10.1007/s10648-013-9244-0
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Aloe14)
#'
#' ## Random-effects meta-analysis
#' meta1 <- meta(cbind(EE,DP,PA),
#'               cbind(V_EE, C_EE_DP, C_EE_PA, V_DP, C_DP_PA, V_PA),
#'               data=Aloe14)
#' ## Remove error code
#' meta1 <- rerun(meta1)
#'
#' summary(meta1)
#'
#' ## Extract the coefficients for the variance component of the random effects
#' coef1 <- coef(meta1, select="random")
#'
#' ## Convert it into a symmetric matrix by row major
#' my.cov <- vec2symMat(coef1, byrow=TRUE)
#'
#' ## Convert it into a correlation matrix
#' cov2cor(my.cov)
#'
#' ## Plot the multivariate effect sizes
#' plot(meta1)
#' }
#'
NULL





#' Compare Nested Models with Likelihood Ratio Statistic
#'
#' It compares nested models with the likelihood ratio statistic from various
#' objects. It is a wrapper of \code{\link[OpenMx]{mxCompare}}.
#'
#'
#' @name anova
#' @aliases anova.wls anova.meta anova.meta3LFIML anova.reml anova.osmasem anova.osmasem2 anova.mxsem
#' @param object An object or a list of objects of various classes. It will be
#' passed to the \code{base} argument in \code{\link[OpenMx]{mxCompare}}.
#' @param \dots An object or a list of objects of various classes. It will be
#' passed to the \code{comparison} argument in \code{\link[OpenMx]{mxCompare}}.
#' @param all A Boolean value on whether to compare all bases with all
#' comparisons. It will be passed to the \code{all} argument in
#' \code{\link[OpenMx]{mxCompare}}.
#' @return A table of comparisons between the models in base and comparison.
#' @note When the objects are class \code{\link[metaSEM]{wls}}, the degrees of
#' freedom in the base and comparison models are incorrect, while the degrees
#' of freedom of the difference between them is correct. If users want to
#' obtain the correct degrees of freedom in the base and comparison models,
#' they may individually apply the \code{\link[metaSEM]{summary}} function on
#' the base and comparison models.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @keywords methods
#' @examples
#'
#' ## Test the significance of a predictor with likelihood ratio test
#' ## Model0: No predictor
#' model0 <- meta(y=yi, v=vi, data=Hox02, model.name="No predictor")
#'
#' ## Model1: With a predictor
#' model1 <- meta(y=yi, v=vi, x=weeks, data=Hox02, model.name="One predictor")
#'
#' ## Compare these two models
#' anova(model1, model0) 
#'
NULL





#' Dataset on the Effectiveness of the BCG Vaccine for Preventing Tuberculosis
#'
#' This dataset includes 13 studies on the effectiveness of the Bacillus
#' Calmette-Guerin (BCG) vaccine for preventing tuberculosis (see van
#' Houwelingen, Arends, & Stijnen (2002) for details).
#'
#' A list of data with the following structure: \describe{ \item{Trial}{Number
#' of the trials} \item{Author}{Authors of the original studies}
#' \item{Year}{Year of publication} \item{VD}{Vaccinated group with disease}
#' \item{VWD}{Vaccinated group without the disease} \item{NVD}{Not vaccinated
#' group with disease} \item{NVWD}{Not vaccinated group without the disease}
#' \item{Latitude}{Geographic latitude of the place where the study was done}
#' \item{Allocation}{Method of treatment allocation} \item{ln_OR}{Natural
#' logarithm of the odds ratio: log((VD/VWD)/(NVD/NVWD))}
#' \item{v_ln_OR}{Sampling variance of ln_OR: 1/VD+1/VWD+1/NVD+1/NVWD}
#' \item{ln_Odd_V}{Natural logarithm of the odds of the vaccinated group:
#' log(VD/VWD)} \item{ln_Odd_NV}{Natural logarithm of the odds of the not
#' vaccinated group: log(NVD/NVWD)} \item{v_ln_Odd_V}{Sampling variance of
#' ln_Odd_V: 1/VD+1/VWD} \item{cov_V_NV}{Sampling covariance between ln_Odd_V
#' and ln_Odd_NV: It is always 0} \item{v_ln_Odd_NV}{Sampling variance of
#' ln_Odd_NV: 1/NVD+1/NVWD} }
#'
#' @name BCG
#' @docType data
#' @references Berkey, C. S., Hoaglin, D. C., Mosteller, F., & Colditz, G. A.
#' (1995). A random-effects regression model for meta-analysis.
#' \emph{Statistics in Medicine}, \bold{14}, 395--411.
#'
#' van Houwelingen, H. C., Arends, L. R., & Stijnen, T. (2002). Advanced
#' methods in meta-analysis: Multivariate approach and meta-regression.
#' \emph{Statistics in Medicine}, \bold{21}, 589--624.
#'
#' Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor
#' package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48.
#' \url{https://www.jstatsoft.org/v36/i03/}.
#' @source Colditz, G. A., Brewer, T. F., Berkey, C. S., Wilson, M. E.,
#' Burdick, E., Fineberg, H. V., & Mosteller, F. (1994). Efficacy of BCG
#' vaccine in the prevention of tuberculosis: Meta-analysis of the published
#' literature. \emph{Journal of the American Medical Association}, \bold{271},
#' 698--702.
#' @keywords datasets
#' @examples
#'
#' data(BCG)
#'
#' ## Univariate meta-analysis on the log of the odds ratio
#' summary( meta(y=ln_OR, v=v_ln_OR, data=BCG,
#'               x=cbind(scale(Latitude,scale=FALSE),
#'               scale(Year,scale=FALSE))) )
#'
#' ## Multivariate meta-analysis on the log of the odds
#' ## The conditional sampling covariance is 0
#' bcg <- meta(y=cbind(ln_Odd_V, ln_Odd_NV), data=BCG,
#'             v=cbind(v_ln_Odd_V, cov_V_NV, v_ln_Odd_NV))
#' summary(bcg)
#'
#' plot(bcg)
#'
NULL





#' Ten Studies of Correlation Matrices used by Becker (2009)
#'
#' This dataset includes ten studies on the relationships between CSAI
#' subscales and sports behavior. The original data were used in Craft et al.
#' (2003), whereas a subset of them was illustrated in Becker (2009).
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of 4x4 correlation matrices. The variables are \emph{Performance},
#' \emph{Cognitive}, \emph{Somatic}, and \emph{Self_confidence}} \item{n}{A
#' vector of sample sizes} \item{Type_of_sport}{Samples based on
#' \emph{Individual} or \emph{Team}} }
#'
#' @name Becker09
#' @docType data
#' @references Becker, B. J. (2009). Model-based meta-analysis. In H. Cooper,
#' L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research
#' synthesis and meta-analysis} (2nd ed., pp. 377-395). New York: Russell Sage
#' Foundation.
#' @source Craft, L. L., Magyar, T. M., Becker, B. J., & Feltz, D. L. (2003).
#' The relationship between the Competitive State Anxiety Inventory-2 and sport
#' performance: a meta-analysis. \emph{Journal of Sport and Exercise
#' Psychology}, \bold{25(1)}, 44-65.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Becker09)
#'
#' #### Fixed-effects model
#' ## First stage analysis
#' fixed1 <- tssem1(Becker09$data, Becker09$n, method="FEM")
#' summary(fixed1)
#'
#' ## Model specification in lavaan syntax
#' model <- "## Regression paths
#'           Performance ~ Cog2Per*Cognitive + SO2Per*Somatic + SC2Per*Self_confidence
#'           Self_confidence ~ Cog2SC*Cognitive + SO2SC*Somatic
#'           ## Fix the variances of Cog and SO at 1
#'           Cognitive ~~ 1*Cognitive
#'           Somatic ~~ 1*Somatic
#'           ## Label the correlation between Cog and SO
#'           Cognitive ~~ cor*Somatic
#'           ## Label the error variances of Per and SC
#'           Performance ~~ var_Per*Performance
#'           Self_confidence ~~ var_SC*Self_confidence"
#'
#' ## Display the model
#' plot(model, layout="spring")
#'
#' RAM <- lavaan2RAM(model, obs.variables=c("Performance", "Cognitive",
#'                                          "Somatic", "Self_confidence"))
#' RAM
#'
#' ## ## Equivalent RAM specification using create.mxMatrix()
#' ## A1 <- create.mxMatrix(c(0, "0.1*Cog2Per", "0.1*SO2Per", "0.1*SC2Per",
#' ##                         0, 0, 0, 0,
#' ##                         0, 0, 0, 0,
#' ##                         0, "0.1*Cog2SC", "0.1*SO2SC",0),
#' ##                       type="Full", byrow=TRUE, ncol=4, nrow=4,
#' ##                       as.mxMatrix=FALSE)
#' ## dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- c("Performance", "Cognitive",
#' ##                                             "Somatic", "Self_confidence")
#' ## A1
#' ##
#' ## S1 <- create.mxMatrix(c("0.1*var_Per",
#' ##                         0, 1,
#' ##                         0, "0.1*cor", 1,
#' ##                         0, 0, 0, "0.1*var_SC"), byrow=TRUE, type="Symm",
#' ##                       as.mxMatrix=FALSE)
#' ## dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- c("Performance", "Cognitive",
#' ##                                             "Somatic", "Self_confidence")
#' ## S1
#'
#' ## Second stage analysis
#' fixed2 <- tssem2(fixed1, RAM=RAM, diag.constraints=TRUE,
#'                  intervals.type="LB", model.name="TSSEM2 Becker09",
#'                  mx.algebras=list( Cog=mxAlgebra(Cog2SC*SC2Per, name="Cog"),
#'                                    SO=mxAlgebra(SO2SC*SC2Per, name="SO"),
#'                                    Cog_SO=mxAlgebra(Cog2SC*SC2Per+SO2SC*SC2Per,
#'                                    name="Cog_SO")) )
#' summary(fixed2)
#'
#' ## Display the model with the parameter estimates
#' plot(fixed2, layout="spring")
#'
#' #### Fixed-effects model: with type of sport as cluster
#' ## First stage analysis
#' cluster1 <- tssem1(Becker09$data, Becker09$n, method="FEM",
#'                    cluster=Becker09$Type_of_sport)
#' summary(cluster1)
#'
#' ## Second stage analysis
#' cluster2 <- tssem2(cluster1, RAM=RAM, diag.constraints=TRUE,
#'                  intervals.type="LB", model.name="TSSEM2 Becker09",
#'                  mx.algebras=list( Cog=mxAlgebra(Cog2SC*SC2Per, name="Cog"),
#'                                    SO=mxAlgebra(SO2SC*SC2Per, name="SO"),
#'                                    Cog_SO=mxAlgebra(Cog2SC*SC2Per+SO2SC*SC2Per,
#'                                    name="Cog_SO")) )
#' summary(cluster2)
#'
#' ## Convert the model to semPlotModel object with 2 plots
#' ## Use the short forms of the variable names
#' my.plots <- lapply(X=cluster2, FUN=meta2semPlot, manNames=c("Per","Cog","SO","SC") )
#'
#' ## Load the library
#' library("semPlot")
#'
#' ## Setup two plots
#' layout(t(1:2))
#' ## The labels are overlapped. We may modify it by using layout="spring"
#' semPaths(my.plots[[1]], whatLabels="est", nCharNodes=10, color="orange",
#'          layout="spring", edge.label.cex=0.8)
#' title("Individual sport")
#' semPaths(my.plots[[2]], whatLabels="est", nCharNodes=10, color="skyblue",
#'          layout="spring", edge.label.cex=0.8)
#' title("Team sport")
#'
#'
#' #### Random-effects model
#' ## First stage analysis
#' random1 <- tssem1(Becker09$data, Becker09$n, method="REM", RE.type="Diag")
#' summary(random1)
#'
#' ## Second stage analysis
#' random2 <- tssem2(random1, RAM=RAM, diag.constraints=TRUE,
#'                   intervals.type="LB", model.name="TSSEM2 Becker09",
#'                   mx.algebras=list( Cog=mxAlgebra(Cog2SC*SC2Per, name="Cog"),
#'                                     SO=mxAlgebra(SO2SC*SC2Per, name="SO"),
#'                                     Cog_SO=mxAlgebra(Cog2SC*SC2Per+SO2SC*SC2Per,
#'                                     name="Cog_SO")) )
#' summary(random2)
#'
#' ## Display the model
#' plot(random2, what="path", layout="spring")
#'
#' ## Display the model with the parameter estimates
#' plot(random2, layout="spring", color="yellow")
#'
#' #### Univariate r approach
#' #### First stage of the analysis
#' uni1 <- uniR1(Becker09$data, Becker09$n)
#' uni1
#'
#' #### Second stage of analysis using OpenMx
#' model2 <- "## Regression paths
#'            Performance ~ Cog2Per*Cognitive + SO2Per*Somatic + SC2Per*Self_confidence
#'            Self_confidence ~ Cog2SC*Cognitive + SO2SC*Somatic
#'            ## Provide starting values for Cog and SO
#'            Cognitive ~~ start(1)*Cognitive
#'            Somatic ~~ start(1)*Somatic
#'            ## Label the correlation between Cog and SO
#'            Cognitive ~~ cor*Somatic
#'            ## Label the error variances of Per and SC
#'            Performance ~~ var_Per*Performance
#'            Self_confidence ~~ var_SC*Self_confidence"
#'
#' RAM2 <- lavaan2RAM(model2, obs.variables=c("Performance", "Cognitive",
#'                                            "Somatic", "Self_confidence"))
#' RAM2
#'
#' uni2mx <- uniR2mx(uni1, RAM=RAM2)
#' summary(uni2mx)
#'
#' #### Second stage of analysis Using lavaan
#' model3 <- "## Regression paths
#'            Performance ~ Cognitive + Somatic + Self_confidence
#'            Self_confidence ~ Cognitive + Somatic"
#'
#' uni2lavaan <- uniR2lavaan(uni1, model3)
#' lavaan::summary(uni2lavaan)
#' }
#'
NULL





#' Studies on Sex Differences in Conformity Reported by Becker (1983)
#'
#' The data set includes studies on sex differences in conformity using the
#' fictitious norm group paradigm reported by Becker (1983).
#'
#' The variables are: \describe{ \item{study}{study number}
#' \item{di}{Standardized mean difference} \item{vi}{Sampling variance of the
#' effect size} \item{percentage}{Percentage of male authors}
#' \item{items}{Number of items} }
#'
#' @name Becker83
#' @docType data
#' @references Cheung, M. W.-L. (2010). Fixed-effects meta-analyses as
#' multiple-group structural equation models. \emph{Structural Equation
#' Modeling}, \bold{17}, 481-509.
#' @source Becker, B. J. (1983, April). Influence again: A comparison of
#' methods for meta-analysis. \emph{Paper presented at the annual meeting of
#' the American Educational Research Association, Montreal.}
#'
#' Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for
#' meta-analysis.} Orlando, FL: Academic Press.
#' @keywords datasets
#' @examples
#'
#' data(Becker83)
#'
#' ## Random-effects meta-analysis
#' summary( meta(y=di, v=vi, data=Becker83) )
#'
#' ## Mixed-effects meta-analysis with log(items) as the predictor
#' summary( meta(y=di, v=vi, x=log(items), data=Becker83) )
#'
NULL





#' Six Studies of Correlation Matrices reported by Becker (1992; 1995)
#'
#' This data set includes six studies of correlation matrices reported by
#' Becker (1992; 1995).
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of 6 studies of correlation matrices. The variables are \emph{Math} (math
#' aptitude), \emph{Spatial} (spatial ability), and \emph{Verbal} (verbal
#' ability)} \item{n}{A vector of sample sizes} }
#'
#' @name Becker92
#' @docType data
#' @source Becker, B. J. (1992). Using results from replicated studies to
#' estimate linear models. \emph{Journal of Educational Statistics},
#' \bold{17(4)}, 341-362. doi:10.3102/10769986017004341
#'
#' Becker, B. J. (1995). Corrections to "Using Results from Replicated Studies
#' to Estimate Linear Models." \emph{Journal of Educational and Behavioral
#' Statistics}, \bold{20(1)}, 100-102. doi:10.2307/1165390
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Becker92)
#'
#' #### Fixed-effects model
#' ## First stage analysis
#' ## Replicate Becker's (1992) analysis using 4 studies only
#'
#' fixed1 <- tssem1(Becker92$data[1:4], Becker92$n[1:4], method="FEM")
#' summary(fixed1)
#'
#' ## ## Prepare a regression model using create.mxMatrix()
#' ## A1 <- create.mxMatrix(c(0,0,0,"0.2*Spatial2Math",
#' ##                         0,0,"0.2*Verbal2Math",0,0), type="Full",
#' ##                         ncol=3, nrow=3, as.mxMatrix=FALSE)
#'
#' ## var.names <- c("Math_aptitude","Spatial","Verbal") 
#'
#' ## ## This step is not necessary but it is useful for inspecting the model.
#' ## dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- var.names
#'
#' ## ## Display A1
#' ## A1
#'
#' ## S1 <- create.mxMatrix(c("0.2*ErrorVarMath",0,0,1,"0.2*CorSpatialVerbal",1),
#' ##                         type="Symm", as.mxMatrix=FALSE)
#'
#' ## ## This step is not necessary but it is useful for inspecting the model.
#' ## dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- var.names
#'
#' ## ## Display S1
#' ## S1
#'
#' ################################################################################
#' ## Alternative model specification in lavaan model syntax
#' model <- "## Regression paths
#'           Math ~ Spatial2Math*Spatial + Verbal2Math*Verbal
#'           Spatial ~~ CorSpatialVerbal*Verbal
#'           ## Fix the variances of Spatial and Verbal at 1
#'           Spatial ~~ 1*Spatial
#'           Verbal ~~ 1*Verbal
#'           ## Label the error variance of Math
#'           Math ~~ ErrorVarMath*Math + start(0.2)*Math"
#'
#' ## Display the model
#' plot(model)
#'
#' RAM <- lavaan2RAM(model, obs.variables=c("Math", "Spatial", "Verbal"))
#' RAM
#'
#' ################################################################################
#'
#' ## Fixed-effects model: Second stage analysis
#' ## Two equivalent versions to calculate the R2 and its 95% LBCI
#' fixed2 <- tssem2(fixed1, RAM=RAM, intervals.type="LB",
#'                  mx.algebras=list(R1=mxAlgebra(Spatial2Math^2+Verbal2Math^2
#'                         +2*CorSpatialVerbal*Spatial2Math*Verbal2Math, name="R1"),
#'                         R2=mxAlgebra(One-Smatrix[1,1], name="R2"),
#'                         One=mxMatrix("Iden", ncol=1, nrow=1, name="One")))
#' summary(fixed2)
#'
#' ## Display the model with the parameter estimates
#' plot(fixed2)
#'
#' #### Random-effects model
#' ## First stage analysis
#' ## No random effects for off-diagonal elements
#' random1 <- tssem1(Becker92$data, Becker92$n, method="REM", RE.type="Diag")
#' summary(random1)
#'
#' ## Random-effects model: Second stage analysis
#' random2 <- tssem2(random1, RAM=RAM)
#' summary(random2)
#'
#' ## Display the model with the parameter estimates
#' plot(random2, color="yellow")
#'
#' #### Similar to conventional fixed-effects GLS approach
#' ## First stage analysis
#' ## No random effects
#' ## Replicate Becker's (1992) analysis using 4 studies only
#' gls1 <- tssem1(Becker92$data[1:4], Becker92$n[1:4], method="REM", RE.type="Zero",
#'                model.name="Fixed effects GLS Stage 1")
#' summary(gls1)
#'
#' ## Fixed-effects GLS model: Second stage analysis
#' gls2 <- tssem2(gls1, RAM=RAM, model.name="Fixed effects GLS Stage 2")
#' summary(gls2)
#' }
#'
NULL





#' Five Studies of Ten Correlation Matrices reported by Becker and Schram
#' (1994)
#'
#' This data set includes five studies of ten correlation matrices reported by
#' Becker and Schram (1994).
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of 10 correlation matrices. The variables are \emph{Math} (math aptitude),
#' \emph{Spatial} (spatial ability), and \emph{Verbal} (verbal ability)}
#' \item{n}{A vector of sample sizes} \item{gender}{\emph{Females} or
#' \emph{Males} samples} }
#'
#' @name Becker94
#' @docType data
#' @source Becker, B. J., & Schram, C. M. (1994). Examining explanatory models
#' through research synthesis. In H. Cooper & L. V. Hedges (Eds.), \emph{The
#' handbook of research synthesis} (pp. 357-381). New York: Russell Sage
#' Foundation.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Becker94)
#'
#' #### Fixed-effects model
#' ## First stage analysis
#' fixed1 <- tssem1(Becker94$data, Becker94$n, method="FEM")
#' summary(fixed1)
#'
#' ## Prepare a regression model using create.mxMatrix()
#' ## A1 <- create.mxMatrix(c(0,0,0,"0.2*Spatial2Math",
#' ##                         0,0,"0.2*Verbal2Math",0,0), type="Full",
#' ##                       ncol=3, nrow=3, name="A1")
#' ## S1 <- create.mxMatrix(c("0.2*ErrorVarMath",0,0,1,
#' ##                         "0.2*CorBetweenSpatialVerbal",1),
#' ##                       type="Symm", name="S1")
#'
#' ## An alternative method to create a regression model with the lavaan syntax
#' model <- "## Regression model
#'           Math ~ Spatial2Math*Spatial + Verbal2Math*Verbal
#'           ## Error variance of Math
#'           Math ~~ ErrorVarMath*Math
#'           ## Variances of Spatial and Verbal fixed at 1.0
#'           Spatial ~~ 1*Spatial
#'           Verbal ~~ 1*Verbal
#'           ## Correlation between Spatial and Verbal
#'           Spatial ~~ CorBetweenSpatialVerbal*Verbal"
#'
#' ## Display the model
#' plot(model)
#'
#' RAM <- lavaan2RAM(model, obs.variables=c("Math", "Spatial", "Verbal"))
#' RAM
#'
#' ## Second stage analysis
#' ## A1 <- RAM$A
#' ## S1 <- RAM$S
#' ## fixed2 <- tssem2(fixed1, Amatrix=A1, Smatrix=S1, intervals.type="LB")
#'
#' fixed2 <- tssem2(fixed1, RAM=RAM, intervals.type="LB")
#' summary(fixed2)
#'
#' ## Display the model with the parameter estimates
#' plot(fixed2)
#'
#' #### Fixed-effects model: with gender as cluster
#' ## First stage analysis
#' cluster1 <- tssem1(Becker94$data, Becker94$n, method="FEM", cluster=Becker94$gender)
#' summary(cluster1)
#'
#' ## Second stage analysis  
#' cluster2 <- tssem2(cluster1, RAM=RAM, intervals.type="LB")
#' summary(cluster2)
#'
#' #### Conventional fixed-effects GLS approach
#' ## First stage analysis
#' ## No random effects
#' ## Replicate Becker's (1992) analysis using 4 studies only
#' gls1 <- tssem1(Becker92$data[1:4], Becker92$n[1:4], method="REM", RE.type="Zero",
#'                model.name="Fixed effects GLS Stage 1")
#' summary(gls1)
#'
#' ## Fixed-effects GLS model: Second stage analysis
#' gls2 <- tssem2(gls1, RAM=RAM, intervals.type="LB",
#'                model.name="Fixed effects GLS Stage 2")
#' summary(gls2)
#' }
#'
NULL





#' Five Published Trails from Berkey et al. (1998)
#'
#' The data set includes five published trials, reported by Berkey et al.
#' (1998), comparing surgical and non-surgical treatments for medium-severity
#' periodontal disease, one year after treatment.
#'
#' The variables are: \describe{ \item{trial}{Trial number}
#' \item{pub_year}{Publication year} \item{no_of_patients}{Number of patients}
#' \item{PD}{Patient improvements (mm) in \emph{probing depth}}
#' \item{AL}{Patient improvements (mm) in \emph{attachment level}}
#' \item{var_PD}{Sampling variance of PD} \item{cov_PD_AL}{Sampling covariance
#' between PD and AD} \item{var_AL}{Sampling variance of AL} }
#'
#' @name Berkey98
#' @docType data
#' @source Berkey, C. S., Hoaglin, D. C., Antczak-Bouckoms, A., Mosteller, F, &
#' Colditz, G. A. (1998). Meta-analysis of multiple outcomes by regression with
#' random effects. \emph{Statistics in Medicine}, \bold{17}, 2537-2550.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Berkey98)
#'
#' #### ML estimation method
#' ## Multivariate meta-analysis
#' x <- meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL), data=Berkey98)
#' x <- rerun(x)
#' summary(x)
#' plot(x)
#'
#' ## Plot individual studies proportional to the weights
#' plot(x, study.weight.plot=TRUE)
#'
#' ## Include forest plot from the metafor package
#' library(metafor)
#' plot(x, diag.panel=TRUE, main="Multivariate meta-analysis",
#'      axis.label=c("PD", "AL"))
#'      forest( rma(yi=PD, vi=var_PD, data=Berkey98) )
#'      title("Forest plot of PD")
#'      forest( rma(yi=AL, vi=var_AL, data=Berkey98) )
#'      title("Forest plot of AL")
#'
#' ## Multivariate meta-analysis with "publication year-1979" as the predictor
#' summary( meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
#'               x=scale(pub_year, center=1979), data=Berkey98,
#'               RE.lbound=NA) )
#'
#' ## Multivariate meta-analysis with equality constraint on the regression coefficients
#' summary( meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
#'               x=scale(pub_year, center=1979), data=Berkey98,
#'               coef.constraints=matrix(c("0.3*Eq_slope", "0.3*Eq_slope"),
#'               nrow=2)) )
#'
#' #### REML estimation method
#' ## Multivariate meta-analysis
#' summary( reml(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
#'               data=Berkey98,
#'               model.name="Multivariate meta analysis with REML") )
#'
#' ## Multivariate meta-analysis with "publication year-1979" as the predictor
#' ## Diagonal structure for the variance component
#' summary( reml(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
#'               RE.constraints=Diag(c("1e-5*Tau2_1_1", "1e-5*Tau2_2_2")),
#'               x=scale(pub_year, center=1979), data=Berkey98) )
#' }
#'
NULL





#' Correlation Matrices from Boer et al. (2016)
#'
#' The data set includes correlation matrices of leader-member exchange in
#' transformational leadership reported by Boer et al. (2016).
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of correlation matrices. The variables are \emph{LMX} (leader-member
#' exchange), \emph{TFL} (transformational leadership), \emph{JS} (job
#' satisfaction), \emph{OC} (organizational commitment), and \emph{LE} (leader
#' effectiveness)} \item{n}{A vector of sample sizes} \item{RelLMX}{The
#' reliability of \emph{LMX}} \item{RelTFL}{The reliability of \emph{TFL}} }
#'
#' @name Boer16
#' @docType data
#' @source Boer, D., Deinert, A., Homan, A. C., & Voelpel, S. C. (2016).
#' Revisiting the mediating role of leader-member exchange in transformational
#' leadership: the differential impact model. \emph{European Journal of Work
#' and Organizational Psychology}, \bold{25}(6), 883-899.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' ## Stage 1 analysis
#' rand1 <- tssem1(Boer16$data, Boer16$n, method="REM", RE.type="Diag",
#'                 acov="weighted")
#' summary(rand1)
#'
#' ## Stage 2 analysis
#' model2a <- 'JS+OC+LE ~ LMX+TFL
#'             LMX ~ TFL
#'             ## Variance of TFL is fixed at 1
#'             TFL ~~ 1*TFL
#'             ## Correlated residuals
#'             JS ~~ OC
#'             JS ~~ LE
#'             OC ~~ LE'
#'
#' ## Display the model
#' plot(model2a)    
#'
#' RAM2a <- lavaan2RAM(model2a, obs.variables = c("LMX", "TFL", "JS", "OC", "LE"),
#'                     A.notation="on", S.notation="with")
#'
#' rand2a <- tssem2(rand1, RAM=RAM2a)
#' summary(rand2a)
#'
#' ## Display the model with the parameter estimates
#' plot(rand2a, layout="spring")    
#' }
#'
NULL





#' A Dataset from Bornmann et al. (2007)
#'
#' A dataset from Bornmann et al. (2007) for three-level meta-analysis.
#'
#' The variables are: \describe{ \item{ID}{ID of the study} \item{Study}{Study
#' name} \item{Cluster}{Cluster for effect sizes} \item{logOR}{Effect size: log
#' odds ratio} \item{v}{Sampling variance of logOR} \item{Year}{Year of
#' publication} \item{Type}{Type of proposal: either \bold{Grant} or
#' \bold{Fellowship}} \item{Discipline}{Discipline of the proposal: either
#' \bold{Physical sciences}, \bold{Life sciences/biology}, \bold{Social
#' sciences/humanities} or \bold{Multidisciplinary})} \item{Country}{Country of
#' the proposal: either the \bold{United States}, \bold{Canada},
#' \bold{Australia}, \bold{United Kingdom} or \bold{Europe}} }
#'
#' @name Bornmann07
#' @docType data
#' @references Cheung, M. W.-L. (2014). Modeling dependent effect sizes with
#' three-level meta-analyses: A structural equation modeling approach.
#' \emph{Psychological Methods}, \bold{19}, 211-229.
#'
#' Marsh, H. W., Bornmann, L., Mutz, R., Daniel, H.-D., & O'Mara, A. (2009).
#' Gender Effects in the Peer Reviews of Grant Proposals: A Comprehensive
#' Meta-Analysis Comparing Traditional and Multilevel Approaches. \emph{Review
#' of Educational Research}, \bold{79(3)}, 1290-1326.
#' doi:10.3102/0034654309334143
#' @source Bornmann, L., Mutz, R., & Daniel, H.-D. (2007). Gender differences
#' in grant peer review: A meta-analysis. \emph{Journal of Informetrics},
#' \bold{1(3)}, 226-238. doi:10.1016/j.joi.2007.03.001
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Bornmann07)
#'
#' #### ML estimation method
#' ## No predictor
#' summary( meta3L(y=logOR, v=v, cluster=Cluster, data=Bornmann07) )
#'
#' ## Type as a predictor
#' ## Grant: 0
#' ## Fellowship: 1
#' summary( meta3L(y=logOR, v=v, x=(as.numeric(Type)-1),
#'                cluster=Cluster, data=Bornmann07) )
#'
#' ## Centered Year as a predictor
#' summary( meta3L(y=logOR, v=v, x=scale(Year, scale=FALSE),
#'                cluster=Cluster, data=Bornmann07) )
#'
#' #### REML estimation method
#' ## No predictor
#' summary( reml3L(y=logOR, v=v, cluster=Cluster, data=Bornmann07) )
#'
#' ## Type as a predictor
#' ## Grants: 0
#' ## Fellowship: 1
#' summary( reml3L(y=logOR, v=v, x=(as.numeric(Type)-1),
#'                 cluster=Cluster, data=Bornmann07) )
#'
#' ## Centered Year as a predictor
#' summary( reml3L(y=logOR, v=v, x=scale(Year, scale=FALSE),
#'                 cluster=Cluster, data=Bornmann07) )
#'
#' ## Handling missing covariates with FIML
#' ## MCAR
#' ## Set seed for replication
#' set.seed(1000000)
#'
#' ## Copy Bornmann07 to my.df
#' my.df <- Bornmann07
#' ## "Fellowship": 1; "Grant": 0
#' my.df$Type_MCAR <- ifelse(Bornmann07$Type=="Fellowship", yes=1, no=0)
#'
#' ## Create 17 out of 66 missingness with MCAR
#' my.df$Type_MCAR[sample(1:66, 17)] <- NA
#' summary(meta3LFIML(y=logOR, v=v, cluster=Cluster, x2=Type_MCAR, data=my.df))
#'
#' ## MAR
#' Type_MAR <- ifelse(Bornmann07$Type=="Fellowship", yes=1, no=0)
#'
#' ## Create 27 out of 66 missingness with MAR for cases Year<1996
#' index_MAR <- ifelse(Bornmann07$Year<1996, yes=TRUE, no=FALSE)
#' Type_MAR[index_MAR] <- NA
#'
#' ## Include auxiliary variable
#' summary(meta3LFIML(y=logOR, v=v, cluster=Cluster, x2=Type_MAR, av2=Year, data=my.df))
#' }
#'
NULL





#' Dataset from Chan, Jones, Jamieson, and Albarracin (2017)
#'
#' A dataset of multiple treatment effects of standardized mean differences on
#' misinformation and debunking effects.
#'
#' The sampling variances and covariances are calculated using Gleser and
#' Olkin's (2009) method for multiple treatment effects (Equations 3.3 and
#' 3.4). Since the sample sizes of the misinformation, debunking, and control
#' groups are not given, it is assumed they are equal.
#'
#' @name Chan17
#' @docType data
#' @format A data frame with 34 independent samples from 6 research reports.
#' \describe{ \item{\code{Author}}{a character vector of study}
#' \item{\code{g_misinfo}}{Hedges' g of misinformation comparing the
#' misinformation experimental and control groups}
#' \item{\code{g_debunk}}{Hedges' g of debunking comparing the debuking
#' experimental and misinformation experimental groups}
#' \item{\code{v_misinfo}}{sampling variance of g_misinfo}
#' \item{\code{c_mis_deb}}{Sampling covariance between \code{g_misinfo} and
#' \code{g_debunk} due to the overlap of the misinformation experimental group}
#' \item{\code{v_debunk}}{sampling variance of g_debunk}
#' \item{\code{PublicationYear}}{publication year}
#' \item{\code{Published}}{published or unpublished}
#' \item{\code{MeanAge}}{mean age of participants}
#' \item{\code{PctFemale}}{percentage of female participants} }
#' @references Gleser, L. J., & Olkin, I. (2009). Stochastically dependent
#' effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.),
#' \emph{The handbook of research synthesis and meta-analysis.} (2nd ed., pp.
#' 357-376). Russell Sage Foundation.
#' @source Chan, M. S., Jones, C. R., Hall Jamieson, K., & Albarracin, D.
#' (2017). Debunking: A meta-analysis of the psychological efficacy of messages
#' countering misinformation. \emph{Psychological Science}, \bold{28(11)},
#' 1531-1546. https://doi.org/10.1177/0956797617714579
#' @keywords datasets
NULL





#' Fifty Studies of Correlation Matrices used in Cheung and Chan (2000)
#'
#' This data set includes fifty studies of correlation matrices on the theory
#' of planned theory reported by Cheung and Chan (2000).
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of 50 studies of correlation matrices. The variables are the attitude toward
#' behavior \emph{att}, subjective norm \emph{sn}, behavioral intention
#' \emph{bi}, and behavior \emph{beh}} \item{n}{A vector of sample sizes} }
#'
#' @name Cheung00
#' @docType data
#' @note These studies were extracted from the original data set for
#' illustration purpose. Some samples contained two or more correlation
#' matrices, and only one of them was arbitrarily selected to avoid the problem
#' of dependence. Moreover, studies with less than 3 correlation coefficients
#' were also excluded.
#' @references Cheung, M.W.-L., & Cheung, S.-F. (2016). Random-effects models
#' for meta-analytic structural equation modeling: Review, issues, and
#' illustrations. \emph{Research Synthesis Methods}, \bold{7}, 140-155.
#' @source Cheung, S.-F., & Chan, D. K.-S. (2000). The role of perceived
#' behavioral control in predicting human behavior: A meta-analytic review of
#' studies on the theory of planned behavior. \emph{Unpublished manuscript},
#' Chinese University of Hong Kong.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Cheung00)
#'
#' ## Full mediation model in lavaan syntax
#' model <- "## Regression paths
#'           bi ~ att2bi*att + sn2bi*sn
#'           beh ~ bi2beh*bi
#'           ## Variances of att and sn are fixed at 1
#'           att ~~ 1*att
#'           sn ~~ 1*sn
#'           ## Covariance between att and sn
#'           att ~~ cov_att_sn*sn
#'           ## Error variances
#'           bi ~~ e_bi*bi
#'           beh ~~ e_beh*beh"
#'
#' RAM <- lavaan2RAM(model, obs.variables=colnames(Cheung00$data[[1]]))
#' RAM
#'
#' ## ## Equivalent RAM specification
#' ## labels <- colnames(Cheung00$data[[1]])
#' ## A <- matrix(c("0","0","0","0",
#' ##               "0","0","0","0",
#' ##               ".2*att2bi", ".2*sn2bi", "0", "0",
#' ##               "0", "0", ".2*bi2beh", "0"),
#' ##             byrow=TRUE, 4, 4)
#' ## dimnames(A) <- list(labels, labels)
#' ## A
#' ##
#' ## S <- create.mxMatrix(c("1",
#' ##                        ".2*cov_att_sn", "1",
#' ##                        0, 0, ".2*e_bi",
#' ##                        0, 0, 0, ".2*e_beh"),
#' ##                      type="Symm", as.mxMatrix=FALSE, byrow=TRUE)
#' ## dimnames(S) <- list(labels, labels)
#' ## S
#'
#' #### Random-effects model
#'
#' ## Stage 1 analysis
#' random_1 <- tssem1(Cheung00$data, Cheung00$n, method="REM", RE.type="Diag",
#'                    acov="weighted")
#' summary(random_1)
#'
#' ## Stage 2 analysis
#' random_2 <- tssem2(random_1, RAM=RAM, intervals.type="LB",
#'                    diag.constraints=TRUE)
#' summary(random_2)
#'
#' ## Display the model
#' plot(random_2, what="path")
#'
#' ## Display the model with the parameter estimates
#' plot(random_2, color="yellow")
#'
#' ## Load the library
#' library("semPlot")
#' }
#'
NULL





#' A Dataset from TSSEM User's Guide Version 1.11 by Cheung (2009)
#'
#' Four studies were selected from the data set used by Cheung and Chan (2005;
#' 2009). Some variables were randomly deleted to illustrate the analysis with
#' missing data.
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of 4 studies of correlation matrices } \item{n}{A vector of sample sizes} }
#'
#' @name Cheung09
#' @docType data
#' @references Cheung, M. W.-L., & Chan, W. (2005). Meta-analytic structural
#' equation modeling: A two-stage approach. \emph{Psychological Methods},
#' \bold{10}, 40-64.
#'
#' Cheung, M. W.-L., & Chan, W. (2009). A two-stage approach to synthesizing
#' covariance matrices in meta-analytic structural equation modeling.
#' \emph{Structural Equation Modeling}, \bold{16}, 28-53.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Cheung09)
#'
#' #### Fixed-effects model: Stage 1 analysis
#' fixed1 <- tssem1(Cheung09$data, Cheung09$n, method="FEM")
#' summary(fixed1)
#'
#' ## Three-factor CFA model in lavaan syntax
#' model <- "## Factor loadings
#'           f1 =~ f1x1*x1 + f1x2*x2 + f1x3*x3
#'           f2 =~ f2x4*x4 + f2x5*x5 + f2x6*x6 + f2x7*x7
#'           f3 =~ f3x8*x8 + f3x9*x9
#'           ## Factor correlations
#'           f1 ~~ corf2f1*f2
#'           f1 ~~ corf3f1*f3
#'           f2 ~~ corf3f2*f3"
#'
#' RAM1 <- lavaan2RAM(model, obs.variables=paste0("x", 1:9))
#' RAM1
#'
#' ## ## Equivalent RAM specification using create.mxMatrix()
#' ## Phi <- create.mxMatrix( c("0.3*corf2f1","0.3*corf3f1","0.3*corf3f2"),
#' ##                         type="Stand", as.mxMatrix=FALSE )
#' ## Psi <- create.mxMatrix( paste("0.2*e", 1:9, sep=""), type="Diag",
#' ##                         as.mxMatrix=FALSE )
#' ## S1 <- bdiagMat(list(Psi, Phi))
#' ## S1 <- as.mxMatrix(S1)
#' ##
#' ## Lambda <- create.mxMatrix( c(".3*f1x1",".3*f1x2",".3*f1x3",rep(0,9),
#' ##                              ".3*f2x4",".3*f2x5",".3*f2x6",".3*f2x7",
#' ##                              rep(0,9),".3*f3x8",".3*f3x9"), type="Full",
#' ##                              ncol=3, nrow=9, as.mxMatrix=FALSE )
#' ## Zero1 <- matrix(0, nrow=9, ncol=9)
#' ## Zero2 <- matrix(0, nrow=3, ncol=12)
#' ## A1 <- rbind( cbind(Zero1, Lambda), Zero2 )
#' ## A1 <- as.mxMatrix(A1)
#' ##
#' ## F1 <- create.Fmatrix(c(rep(1,9), rep(0,3)))
#'
#' #### Fixed-effects model: Stage 2 analysis
#' fixed2 <- tssem2(fixed1, RAM=RAM1, intervals.type="LB")
#' summary(fixed2)
#'
#' ## Display the model
#' plot(fixed2, what="path")
#'
#' ## Display the model with the parameter estimates
#' plot(fixed2, latNames=c("f1", "f2", "f3"), edge.label.cex=0.8,
#'      color="yellow")
#' }
#'
NULL





#' Extract Parameter Estimates from various classes.
#'
#' It extracts the parameter estimates from objects of various classes.
#'
#'
#' @name coef
#' @aliases coef.tssem1FEM coef.tssem1FEM.cluster coef.tssem1REM coef.wls coef.wls.cluster coef.meta coef.meta3LFIML coef.reml coef.osmasem coef.osmasem2 coef.mxsem
#' @param object An object returned from either class \code{tssem1FEM}, class
#' \code{tssem1FEM.cluster}, class \code{tssem1REM}, class \code{wls}, class
#' \code{wls.cluster}, class \code{meta}, class \code{reml}, class
#' \code{osmasem}, class \code{osmasem2}, or class \code{sem}
#' @param select Select \code{all} for both fixed- and random-effects
#' parameters, \code{fixed} for the fixed-effects parameters or \code{random}
#' for the random-effects parameters. For \code{meta3LFIML} objects,
#' \code{allX} is used to extract all parameters including the predictors and
#' auxiliary variables.
#' @param \dots Further arguments; currently none is used
#' @return Parameter estimates for both fixed-effects (if any) and
#' random-effects (if any)
#' @note \code{coef.sem} is simply a wraper of \code{omxGetParameters}. Extra
#' arguments will be passed to it
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
#' \code{\link[metaSEM]{meta}}, \code{\link[metaSEM]{reml}},
#' \code{\link[OpenMx]{omxGetParameters}}, \code{\link[metaSEM]{osmasem}}
#' @keywords methods
#' @examples
#'
#' ## Random-effects meta-analysis
#' model1 <- meta(y=yi, v=vi, data=Hox02)
#' coef(model1)
#'
#' ## Fixed-effects only
#' coef(model1, select="fixed")
#'
NULL





#' Correlation Matrices from Cooke et al. (2016)
#'
#' The data set includes correlation matrices on using the theory of planned
#' behavior to predict alcohol consumption reported by Cooke et al. (2016).
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of correlation matrices. The variables are \emph{SN} (subjective norm),
#' \emph{ATT} (attitude), \emph{PBC} (perceived behavior control), \emph{BI}
#' (behavioral intention), and \emph{BEH} (behavior).} \item{n}{A vector of
#' sample sizes.} \item{MeanAge}{Mean age of the participants except for
#' \code{Ajzen and Sheikh (2013)}, which is the median age, and \code{Glassman,
#' et al. (2010a)} to \code{Glassman, et al. (2010d)}, which are based on the
#' range of 18 to 24.} \item{Female}{Percentage of female participants.} }
#'
#' @name Cooke16
#' @docType data
#' @references Cheung, M. W.-L., & Hong, R. Y. (2017). Applications of
#' meta-analytic structural equation modeling in health psychology: Examples,
#' issues, and recommendations. \emph{Health Psychology Review}, \bold{11},
#' 265-279.
#' @source Cooke, R., Dahdah, M., Norman, P., & French, D. P. (2016). How well
#' does the theory of planned behaviour predict alcohol consumption? A
#' systematic review and meta-analysis. \emph{Health Psychology Review},
#' \bold{10}(2), 148-167.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' ## Check whether the correlation matrices are valid (positive definite)
#' Cooke16$data[is.pd(Cooke16$data)==FALSE]
#'
#' ## Since the correlation matrix in Study 3 is not positive definite,
#' ## we exclude it in the following analyses
#' my.data <- Cooke16$data[-3]
#' my.n <- Cooke16$n[-3]
#'
#' ## Show the no. of studies per correlation
#' pattern.na(my.data, show.na = FALSE)
#'
#' ## Show the total sample sizes per correlation
#' pattern.n(my.data, my.n)
#'
#' ## Stage 1 analysis
#' ## Random-effects model
#' random1 <- tssem1(my.data, my.n, method="REM", RE.type="Diag", acov="weighted")
#' summary(random1)
#'
#' ## Model specification in lavaan syntax
#' model <- "## Regression paths
#'           BI ~ SN2BI*SN + ATT2BI*ATT + PBC2BI*PBC
#'           BEH ~ PBC2BEH*PBC + BI2BEH*BI
#'           ## Variances of SN, ATT, and PBC are fixed at 1
#'           SN ~~ 1*SN
#'           ATT ~~ 1*ATT
#'           PBC ~~ 1*PBC
#'           ## Covariances among SN, ATT, and PBC
#'           ATT ~~ ATT_SN*SN
#'           PBC ~~ PBC_SN*SN
#'           PBC ~~ PBC_ATT*ATT
#'           ## Error variances of BI and BEH
#'           BI ~~ VarBI*BI
#'           BEH ~~ VarBEH*BEH"
#'
#' RAM1 <- lavaan2RAM(model, obs.variables=colnames(Cooke16$data[[1]]))
#' RAM1
#'
#' ## ## Equivalent RAM specification using create.mxMatrix()
#' ## A1 <- create.mxMatrix(c(0,0,0,0,0,
#' ##                         0,0,0,0,0,
#' ##                         0,0,0,0,0,
#' ##                         "0.2*SN2BI","0.2*ATT2BI","0.2*PBC2BI",0,0,
#' ##                         0,0,"0.2*PBC2BEH","0.2*BI2BEH",0),
#' ##                         type="Full", ncol=5, nrow=5,
#' ##                         byrow=TRUE, as.mxMatrix=FALSE)
#' ## dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- colnames(Cooke16$data[[1]])
#' ## A1
#' ##
#' ## S1 <- create.mxMatrix(c(1,
#' ##                         "0.1*ATT_SN", 1,
#' ##                         "0.1*PBC_SN", "0.1*PBC_ATT", 1,
#' ##                         0, 0, 0, "0.5*VarBI",
#' ##                         0, 0, 0, 0, "0.5*VarBEH"),
#' ##                       type = "Symm", ncol=5, nrow=5,
#' ##                       byrow=TRUE, as.mxMatrix=FALSE)
#' ## dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- colnames(Cooke16$data[[1]])
#' ## S1
#'
#' ## Stage 2 analysis
#' random2 <- tssem2(random1, RAM=RAM1, diag.constraints=FALSE,
#'                   intervals.type="LB")
#' summary(random2)
#'
#' ## Display the model
#' plot(random2, what="path")    
#'
#' ## Display the model with the parameter estimates
#' plot(random2, color="yellow")
#' }
#'
NULL





#' Selected effect sizes from Cooper et al. (2003)
#'
#' Fifty-six effect sizes from 11 districts from Cooper et al. (2003) were
#' reported by Konstantopoulos (2011).
#'
#' The variables are: \describe{ \item{District}{District ID}
#' \item{Study}{Study ID} \item{y}{Effect size} \item{v}{Sampling variance}
#' \item{Year}{Year of publication} }
#'
#' @name Cooper03
#' @docType data
#' @references Konstantopoulos, S. (2011). Fixed effects and variance
#' components estimation in three-level meta-analysis. \emph{Research Synthesis
#' Methods}, \bold{2}, 61-76. doi:10.1002/jrsm.35
#' @source Cooper, H., Valentine, J. C., Charlton, K., & Melson, A. (2003). The
#' Effects of Modified School Calendars on Student Achievement and on School
#' and Community Attitudes. \emph{Review of Educational Research},
#' \bold{73(1)}, 1-52. doi:10.3102/00346543073001001
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Cooper03)
#'
#' #### ML estimation method
#' ## No predictor
#' summary( model1 <- meta3L(y=y, v=v, cluster=District, data=Cooper03) )
#'
#' ## Show all heterogeneity indices and their 95% confidence intervals
#' summary( meta3L(y=y, v=v, cluster=District, data=Cooper03,
#'                intervals.type="LB", I2=c("I2q", "I2hm", "I2am", "ICC")) )
#'
#' ## Year as a predictor
#' summary( meta3L(y=y, v=v, cluster=District, x=scale(Year, scale=FALSE),
#'                data=Cooper03, model.name="Year as a predictor") )
#'
#' ## Equality of level-2 and level-3 heterogeneity
#' summary( model2 <- meta3L(y=y, v=v, cluster=District, data=Cooper03,
#'                          RE2.constraints="0.2*EqTau2",
#'                          RE3.constraints="0.2*EqTau2",
#'                          model.name="Equal Tau2") )
#'
#' ## Compare model2 vs. model1
#' anova(model1, model2)
#'
#' #### REML estimation method
#' ## No predictor
#' summary( reml3L(y=y, v=v, cluster=District, data=Cooper03) )
#'
#' ## Level-2 and level-3 variances are constrained equally 
#' summary( reml3L(y=y, v=v, cluster=District, data=Cooper03,
#'                RE.equal=TRUE, model.name="Equal Tau2") )
#'
#' ## Year as a predictor
#' summary( reml3L(y=y, v=v, cluster=District, x=scale(Year, scale=FALSE),
#'                data=Cooper03, intervals.type="LB") )
#'
#' ## Handling missing covariates with FIML
#' ## Create 20/56 MCAR data in Year
#' set.seed(10000)
#' Year_MCAR <- Cooper03$Year
#' Year_MCAR[sample(56, 20)] <- NA
#' summary( meta3LFIML(y=y, v=v, cluster=District, x2=scale(Year_MCAR, scale=FALSE),
#'                     data=Cooper03, model.name="NA in Year_MCAR") )
#' }
#'
NULL





#' Factor Correlation Matrices of Big Five Model from Digman (1997)
#'
#' The data set includes fourteen studies of the factor correlation matrices of
#' the Five-Factor Model of personality reported by Digman (1997).
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of 14 studies of correlation matrices. The variables are
#' \emph{Agreeableness} (A), \emph{Conscientiousness} (C), \emph{Emotional
#' Stability} (ES), \emph{Extraversion} (E) and \emph{Intellect} (I)}
#' \item{n}{A vector of sample sizes} \item{cluster}{Types of participants of
#' the studies} }
#'
#' @name Digman97
#' @docType data
#' @references Cheung, M. W.-L., & Chan, W. (2005). Classifying correlation
#' matrices into relatively homogeneous subgroups: A cluster analytic approach.
#' \emph{Educational and Psychological Measurement}, \bold{65}, 954-979.
#' @source Digman, J.M. (1997). Higher-order factors of the Big Five.
#' \emph{Journal of Personality and Social Psychology}, \bold{73}, 1246-1256.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' Digman97
#'
#' ##### Fixed-effects TSSEM
#' fixed1 <- tssem1(Digman97$data, Digman97$n, method="FEM")
#' summary(fixed1)
#'
#' ## Two-factor CFA model in lavaan syntax
#' model <- "## Factor loadings
#'           Alpha=~A+C+ES
#'           Beta=~E+I
#'           ## Factor correlation
#'           Alpha~~Beta"
#'
#' ## Display the model
#' plot(model)
#'
#' RAM <- lavaan2RAM(model, obs.variables=c("A","C","ES","E","I"),
#'                   A.notation="on", S.notation="with")
#' RAM
#'
#' ## ## Equivalent RAM specification
#' ## Phi <- matrix(c(1,"0.3*cor","0.3*cor",1), ncol=2, nrow=2)
#' ## Psi <- Diag(c("0.2*e1","0.2*e2","0.2*e3","0.2*e4","0.2*e5"))
#' ## S1 <- bdiagMat(list(Psi, Phi))
#' ## dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- c("A","C","ES","E","I","Alpha","Beta")
#' ## S1
#' ##
#' ## Lambda <- matrix(c(".3*Alpha_A",".3*Alpha_C",".3*Alpha_ES",rep(0,5),
#' ##                    ".3*Beta_E",".3*Beta_I"), ncol=2, nrow=5)
#' ## A1 <- rbind( cbind(matrix(0,ncol=5,nrow=5), Lambda),
#' ##              matrix(0, ncol=7, nrow=2) )
#' ## dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- c("A","C","ES","E","I","Alpha","Beta")
#' ## A1
#' ##
#' ## F1 <- create.Fmatrix(c(1,1,1,1,1,0,0), as.mxMatrix=FALSE)
#'
#' fixed2 <- tssem2(fixed1, RAM=RAM,
#'                  model.name="TSSEM2 Digman97")
#' summary(fixed2)
#'
#' ## Display the model with the parameter estimates
#' plot(fixed2)
#'
#' #### Fixed-effects TSSEM with several clusters
#' #### Create a variable for different samples
#' #### Younger participants: Children and Adolescents
#' #### Older participants: others
#' cluster <- ifelse(Digman97$cluster %in% c("Children","Adolescents"),
#'                   yes="Younger participants", no="Older participants")
#'
#' #### Show the cluster
#' cluster
#'
#' ## Example of Fixed-effects TSSEM with several clusters
#' fixed1.cluster <- tssem1(Digman97$data, Digman97$n, method="FEM",
#'                          cluster=cluster)
#' summary(fixed1.cluster)
#'
#' fixed2.cluster <- tssem2(fixed1.cluster, RAM=RAM)
#' #### Please note that the estimates for the younger participants are problematic.
#' summary(fixed2.cluster)
#'
#' ## Setup two plots
#' layout(t(1:2))
#'
#' ## Plot the first group
#' plot(fixed2.cluster[[1]])
#' title("Younger participants")
#'
#' ## Plot the second group
#' plot(fixed2.cluster[[2]])
#' title("Older participants")
#'
#' #### Random-effects TSSEM with random effects on the diagonals
#' random1 <- tssem1(Digman97$data, Digman97$n, method="REM",
#'                   RE.type="Diag")
#' summary(random1)
#'
#' random2 <- tssem2(random1, RAM=RAM)
#' summary(random2)
#'
#' ## Display the model with the parameter estimates
#' plot(random2, color="green")
#' }
#'
NULL





#' Two Datasets from Gleser and Olkin (1994)
#'
#' It includes two datasets in multiple-treatment studies and multiple-endpoint
#' studies reported by Gleser and Olkin (1994).
#'
#'
#' @name Gleser94
#' @docType data
#' @format A list of two data frames.  \describe{ \item{\code{MTS}}{A data
#' frame of multiple-treatment studies.} \item{\code{MES}}{A data frame of
#' multiple-endpoint studies.} }
#' @seealso \code{\link[metaSEM]{smdMTS}}, \code{\link[metaSEM]{smdMES}}
#' @source Gleser, L. J., & Olkin, I. (1994). Stochastically dependent effect
#' sizes. In H. Cooper & L. V. Hedges (Eds.), The handbook of research
#' synthesis. (pp. 339-355). New York: Russell Sage Foundation.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Gleser94)
#'
#' #### Multiple-treatment studies
#' Gleser94$MTS
#'
#' ## Assuming homogeneity of variances
#' my.MTS <- t(apply(Gleser94$MTS, MARGIN=1,
#'             function(x)
#'             smdMTS(m=x[c("Mean.C", "Mean.E1", "Mean.E2", "Mean.E3", "Mean.E4", "Mean.E5")],
#'                    v=x[c("SD.C", "SD.E1", "SD.E2", "SD.E3", "SD.E4", "SD.E5")]^2,
#'                    n=x[c("N.C", "N.E1", "N.E2", "N.E3", "N.E4", "N.E5")],
#'                    homogeneity="variance", list.output=FALSE)))
#'
#' ## Fixed-effects multivariate meta-analysis
#' fit.MTS <- meta(y=my.MTS[, 1:5], 
#'                 v=my.MTS[, 6:20], 
#'                 RE.constraints = diag(0, ncol=5, nrow=5),
#'                 model.name="MTS")
#' summary(fit.MTS)
#'
#' #### Multiple-endpoint studies
#' Gleser94$MES
#'
#' ## Calculate the sampling variances and covariance and amend into the data set
#' Gleser94$MES$Uncoached.V11 <- with(Gleser94$MES, SD.Uncoached.Math^2)
#' Gleser94$MES$Uncoached.V21 <- with(Gleser94$MES,
#'                                    SD.Uncoached.Math*Cor.Math.Verbal*SD.Uncoached.Verbal)
#' Gleser94$MES$Uncoached.V22 <- with(Gleser94$MES, SD.Uncoached.Verbal^2)
#'
#' Gleser94$MES$Coached.V11 <- with(Gleser94$MES, SD.Coached.Math^2)
#' Gleser94$MES$Coached.V21 <- with(Gleser94$MES,
#'                                  SD.Coached.Math*Cor.Math.Verbal*SD.Coached.Verbal)
#' Gleser94$MES$Coached.V22 <- with(Gleser94$MES, SD.Coached.Verbal^2)
#'
#' ## Assuming homogeneity of covariance matrices
#' my.MES <- t(apply(Gleser94$MES, MARGIN=1,
#'             function(x)
#'             smdMES(m1=x[c("Mean.Uncoached.Math", "Mean.Uncoached.Verbal")],
#'                    m2=x[c("Mean.Coached.Math", "Mean.Coached.Verbal")],
#'                    V1=vec2symMat(x[c("Uncoached.V11", "Uncoached.V21", "Uncoached.V22")]),
#'                    V2=vec2symMat(x[c("Coached.V11", "Coached.V21", "Coached.V22")]),
#'                    n1=x["N.Uncoached"],
#'                    n2=x["N.Coached"],
#'                    homogeneity="covariance", list.output=FALSE)))
#'
#' ## Fixed-effects multivariate meta-analysis
#' fit.MES <- meta(y=my.MES[, 1:2], 
#'                 v=my.MES[, 3:5], 
#'                 RE.constraints = diag(0, ncol=2, nrow=2),
#'                 model.name="MES")
#' summary(fit.MES)
#' }
#'
NULL





#' Correlation Matrices from Gnambs, Scharl, and Schroeders (2018)
#'
#' The data set includes 113 correlation matrices on the Rosenberg Self-Esteem
#' Scale reported by Gnambs, Scharl, and Schroeders (2018). Thirty-six studies
#' were based on the reported correlation matrices (\code{CorMat=1}) whereas
#' the correlation matrices of the other 77 studies were calculated from the
#' reported factor loadings.
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of 113 correlation matrices. The variable names are from \emph{I1} to
#' \emph{I10}.} \item{n}{A vector of sample sizes.} \item{Year}{The year of
#' publications.} \item{Country}{The country of studies conducted.}
#' \item{Language}{The language used in the studies.}
#' \item{Publication}{Whether the studies were published (1) or unpublished
#' (0).} \item{MeanAge}{Mean age of the participants.}
#' \item{FemaleProp}{Proportion of the female participants.}
#' \item{Individualism}{Individualism score of the country.}
#' \item{CorMat}{Whether the correlation matrices are obtained from the
#' original studies (1) or reproduced from the factor loadings (0).} }
#'
#' @name Gnambs18
#' @docType data
#' @source Gnambs, T., Scharl, A., & Schroeders, U. (2018). The structure of
#' the Rosenberg Self-Esteem Scale. \emph{Zeitschrift Fur Psychologie},
#' \bold{226}(1), 14-29. https://doi.org/10.1027/2151-2604/a000317
#' @keywords datasets
NULL





#' Effects of Open Education Reported by Hedges and Olkin (1985)
#'
#' Effects of open education on attitude toward school and on reading
#' achievement reported by Hedges and Olkin (1985).
#'
#' The variables are: \describe{ \item{study}{Study number}
#' \item{d_att}{Standardized mean difference on \emph{attitude}}
#' \item{d_ach}{Standardized mean difference on \emph{achievement}}
#' \item{var_att}{Sampling variance of the effect size of \emph{attitude}}
#' \item{cov_att_ach}{Sampling covariance between the effect sizes}
#' \item{var_ach}{Sampling variance of the effect size of \emph{achievement}} }
#'
#' @name HedgesOlkin85
#' @docType data
#' @references Cheung, M. W.-L. (2010). Fixed-effects meta-analyses as
#' multiple-group structural equation models. \emph{Structural Equation
#' Modeling}, \bold{17}, 481-509.
#' @source Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for
#' meta-analysis.} Orlando, FL: Academic Press.
#' @keywords datasets
#' @examples
#'
#' data(HedgesOlkin85)
#'
#' ## Fixed-effects meta-analysis
#' summary( meta(y=cbind(d_att, d_ach),
#'               v=cbind(var_att, cov_att_ach, var_ach),
#'               data=HedgesOlkin85,
#'               RE.constraints=matrix(0, nrow=2, ncol=2)) )
#'
NULL





#' Simulated Effect Sizes Reported by Hox (2002)
#'
#' Twenty stimulated studies on standardized mean difference and one continuous
#' study characteristic reported by Hox (2002).
#'
#' The variables are: \describe{ \item{study}{Study number} \item{yi}{Effect
#' size (standardized mean difference)} \item{vi}{Sampling variance of the
#' effect size} \item{weeks}{Duration of the experimental intervention in terms
#' of weeks} }
#'
#' @name Hox02
#' @docType data
#' @references Cheung, M. W.-L. (2008). A model for integrating fixed-,
#' random-, and mixed-effects meta-analyses into structural equation modeling.
#' \emph{Psychological Methods}, \bold{13}, 182-202.
#' @source Hox, J. J. (2002). \emph{Multilevel analysis: Techniques and
#' applications.} Mahwah, N.J.: Lawrence Erlbaum Associates.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Hox02)
#'
#' #### ML estimation method
#' ## Random-effects meta-analysis
#' summary( meta(y=yi, v=vi, data=Hox02, I2=c("I2q", "I2hm"), intervals.type="LB") ) 
#'
#' ## Fixed-effects meta-analysis
#' summary( meta(y=yi, v=vi, data=Hox02, RE.constraints=0,
#'               model.name="Fixed effects model") )
#'
#' ## Mixed-effects meta-analysis with "weeks" as a predictor
#' ## Request likelihood-based CI
#' summary( meta(y=yi, v=vi, x=weeks, data=Hox02, intervals.type="LB",
#'               model.name="Mixed effects meta analysis with LB CI") )
#'
#' #### REML estimation method
#' ## Random-effects meta-analysis with REML
#' summary( VarComp <- reml(y=yi, v=vi, data=Hox02) )
#'
#' ## Extract the variance component
#' VarComp_REML <- matrix( coef(VarComp), ncol=1, nrow=1 )
#'
#' ## Meta-analysis by treating the variance component as fixed
#' summary( meta(y=yi, v=vi, data=Hox02, RE.constraints=VarComp_REML) )
#'
#'
#' ## Mixed-effects meta-analysis with "weeks" as a predictor
#' ## Request Wald CI
#' summary( reml(y=yi, v=vi, x=weeks, intervals.type="z",
#'               data=Hox02, model.name="REML with LB CI") )
#' }
#'
NULL





#' Fourteen Studies of Correlation Matrices reported by Hunter (1983)
#'
#' This dataset includes fourteen studies of Correlation Matrices reported by
#' Hunter (1983)
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of 14 studies of correlation matrices. The variables are \emph{Ability},
#' \emph{Job knowledge}, \emph{Work sample} and \emph{Supervisor rating}}
#' \item{n}{A vector of sample sizes} }
#'
#' @name Hunter83
#' @docType data
#' @source Hunter, J. E. (1983). A causal analysis of cognitive ability, job
#' knowledge, job performance, and supervisor ratings. In F. Landy, S. Zedeck,
#' & J. Cleveland (Eds.), \emph{Performance Measurement and Theory} (pp.
#' 257-266). Hillsdale, NJ: Erlbaum.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Hunter83)
#'
#' #### Fixed-effects model
#' ## First stage analysis
#' fixed1 <- tssem1(Hunter83$data, Hunter83$n, method="FEM",
#'                  model.name="TSSEM1 fixed effects model")
#' summary(fixed1)
#'
#' #### Second stage analysis
#' ## Model without direct effect from Ability to Supervisor
#' ## A1 <- create.mxMatrix(c(0,"0.1*A2J","0.1*A2W",0,0,0,"0.1*J2W","0.1*J2S",
#' ##                         0,0,0,"0.1*W2S",0,0,0,0),
#' ##                         type="Full", ncol=4, nrow=4, as.mxMatrix=FALSE)
#'
#' ## ## This step is not necessary but it is useful for inspecting the model.
#' ## dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- c("Ability","Job","Work","Supervisor") 
#' ## A1
#'
#' ## S1 <- create.mxMatrix(c(1,"0.1*Var_e_J", "0.1*Var_e_W", "0.1*Var_e_S"),
#' ##                       type="Diag", as.mxMatrix=FALSE)
#' ## dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- c("Ability","Job","Work","Supervisor") 
#' ## S1
#'
#' ################################################################################
#' ## Model specification in lavaan model syntax
#' ## The "ind" effect can be defined within the syntax
#' model1 <- "## Regression paths
#'            Job_knowledge ~ A2J*Ability
#'            Work_sample ~ A2W*Ability + J2W*Job_knowledge
#'            Supervisor ~ J2S*Job_knowledge + W2S*Work_sample
#'
#'            ## Fix the variance of Ability at 1
#'            Ability ~~ 1*Ability
#'
#'            ## Label the error variances of the dependent variables
#'            Job_knowledge ~~ VarE_J*Job_knowledge
#'            Work_sample ~~ VarE_W*Work_sample
#'            Supervisor ~~ VarE_S*Supervisor
#'
#'            ## Define an indirect effect
#'            ind := A2J*J2S+A2J*J2W*W2S+A2W*W2S"
#'
#' ## Display the model
#' plot(model1, layout="spring", sizeMan=10)
#'
#' RAM1 <- lavaan2RAM(model1, obs.variables=c("Ability","Job_knowledge",
#'                    "Work_sample","Supervisor"))
#' RAM1
#'
#' ################################################################################
#' fixed2 <- tssem2(fixed1, RAM=RAM1, intervals.type="z",
#'                  diag.constraints=FALSE,
#'                  model.name="TSSEM2 fixed effects model")
#' summary(fixed2)
#'
#' ## Display the model with the parameter estimates
#' plot(fixed2, layout="spring")
#'
#' ## Coefficients
#' coef(fixed2)
#'
#' ## VCOV based on parametric bootstrap
#' vcov(fixed2)
#'
#' #### Random-effects model with diagonal elements only
#' ## First stage analysis
#' random1 <- tssem1(Hunter83$data, Hunter83$n, method="REM", RE.type="Diag", 
#'                   acov="weighted", model.name="TSSEM1 random effects model")
#' summary(random1)
#'
#' model2 <- "## Regression paths
#'            Job_knowledge ~ A2J*Ability
#'            Work_sample ~ A2W*Ability + J2W*Job_knowledge
#'            Supervisor ~ J2S*Job_knowledge + W2S*Work_sample
#'
#'            ## Fix the variance of Ability at 1
#'            Ability ~~ 1*Ability
#'
#'            ## Label the error variances of the dependent variables
#'            Job_knowledge ~~ VarE_J*Job_knowledge
#'            Work_sample ~~ VarE_W*Work_sample
#'            Supervisor ~~ VarE_S*Supervisor"
#'
#' RAM2 <- lavaan2RAM(model2, obs.variables=c("Ability","Job_knowledge",
#'                    "Work_sample","Supervisor"))
#' RAM2
#'
#' ## Second stage analysis
#' ## Model without direct effect from Ability to Supervisor
#'
#' ## The "ind" effect is defined in tssem2().
#' random2 <- tssem2(random1, RAM=RAM2, intervals.type="LB",
#'                   diag.constraints=FALSE,
#'                   mx.algebras=
#'                       list(ind=mxAlgebra(A2J*J2S+A2J*J2W*W2S+A2W*W2S, name="ind")),
#'                   model.name="TSSEM2 random effects model")
#'
#' summary(random2)
#'
#' ## Display the model with the parameter estimates
#' plot(random2, layout="spring")
#' }
#'
NULL





#' A Dataset from ISSP (2005)
#'
#' Thirty-two covariance matrices on work-related attitudes were extracted from
#' the International Social Survey Programme 2005: Work Orientation III (ISSP,
#' 2005). Seven variables were selected for demonstration purposes. They were
#' grouped into three constructs: \emph{Importance of Job Prospects} measured
#' by job security (JP1), high income (JP2), and opportunity for advancement
#' (JP3); \emph{Importance of Job Autonomy} measured by work independently
#' (JA1) and decide time of work (JA2); and \emph{Importance of Contributions
#' to Society} measured by help other people (CS1) and a job useful to society
#' (CS2).
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of 32 covariance matrices } \item{n}{A vector of sample sizes}
#' \item{means}{A matrix of means} \item{pdi}{Hofstede's Power Distance Index}
#' \item{idv}{Hofstede's Individualism} \item{mas}{Hofstede's Masculinity}
#' \item{uai}{Hofstede's Uncertainty Avoidance Index} \item{ltowvs}{Hofstede's
#' Long- Versus Short-Term Orientation} \item{ivr}{Hofstede's Indulgence Versus
#' Restraint} }
#'
#' @name issp05
#' @docType data
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{issp89}}
#' @source ISSP Research Group (2007): \emph{International Social Survey
#' Programme 2005: Work Orientation III (ISSP 2005)}. GESIS Data Archive,
#' Cologne. ZA4350 Data file Version 1.0.0, doi:10.4232/1.4350
#'
#' https://geerthofstede.com/research-and-vsm/dimension-data-matrix/
#' @keywords datasets
#' @examples
#'
#'
#' \donttest{
#' data(issp05)
#'
#' #### TSSEM random-effects model with covariance matrices
#'
#' ## Stage 1 analysis
#' rand1 <- tssem1(issp05$data, issp05$n, method="REM", cor.analysis=FALSE)
#' summary(rand1)
#'
#' ## Proposed model
#' model1 <- "JP =~ JP1 + JP2 + JP3
#'            JA =~ JA1 + JA2
#'            CS =~ CS1 + CS2
#'            JP ~~ JA + CS
#'            JA ~~ CS"
#'
#' ram1 <- lavaan2RAM(model1, obs.variables=c("JP1", "JP2", "JP3", "JA1", "JA2",
#'                                            "CS1", "CS2"))
#'
#' ## Stage 2 analysis
#' rand2 <- tssem2(rand1, RAM=ram1)
#' summary(rand2)
#'
#' plot(rand2)
#'
#' #### OSMASEM with covariance matrices
#' ## Create a data frame for the OSMASEM
#' df <- Cor2DataFrame(issp05$data, n=issp05$n, Means=issp05$means,
#'                     cor.analysis=FALSE)
#'
#' ## Standardize idv
#' idv <- scale(issp05$idv)
#'
#' ## Replace missing values with mean
#' idv[is.na(idv)] <- mean(idv, na.rm=TRUE)
#' df$data$idv <- idv
#'
#' ## No moderator 
#' fit1 <- osmasem2(model.name="No_moderator", RAM=ram1, data=df,
#'                  cor.analysis=FALSE, mean.analysis=FALSE)
#' summary(fit1, fitIndices = TRUE)
#'
#' ## Proposed model with idv as a moderator
#' model2 <- "JP =~ a*JP1 + b*JP2 + c*JP3
#'            JA =~ d*JA1 + e*JA2
#'            CS =~ f*CS1 + g*CS2
#'            JP ~~ JA + CS
#'            JA ~~ CS
#'            a == a0 + a1*data.idv
#'            b == b0 + b1*data.idv
#'            c == c0 + c1*data.idv
#'            d == d0 + d1*data.idv
#'            e == e0 + e1*data.idv
#'            f == f0 + f1*data.idv
#'            g == g0 + g1*data.idv"
#'
#' ram2 <- lavaan2RAM(model2, obs.variables=c("JP1", "JP2", "JP3", "JA1", "JA2",
#'                                            "CS1", "CS2"))
#'
#' fit2 <- osmasem2(RAM=ram2, data=df, cor.analysis=FALSE, mean.analysis=FALSE,
#'                  replace.constraints = TRUE)
#' summary(fit2)
#'
#' ## Compare fit1 and fit2
#' anova(fit2, fit1)
#' }
#'
NULL





#' A Dataset from Cheung and Chan (2005; 2009)
#'
#' Eleven covariance matrices on work-related attitudes were extracted from the
#' Inter-University Consortium for Political and Social Research (1989). Nine
#' variables were selected by Cheung and Chan (2005; 2009) for demonstration
#' purposes. They were grouped into three constructs: \emph{Job Prospects}
#' measured by job security (JP1), income (JP2), and advancement opportunity
#' (JP3); \emph{Job Nature} measured by interesting job (JN1), independent work
#' (JN2), help other people (JN3), and useful to society (JN4); and \emph{Time
#' Demand} measured by flexible working hours (TD1) and lots of leisure time
#' (TD2).
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of 11 studies of covariance matrices } \item{n}{A vector of sample sizes} }
#'
#' @name issp89
#' @docType data
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{issp05}}
#' @references Cheung, M. W.-L., & Chan, W. (2005). Meta-analytic structural
#' equation modeling: A two-stage approach. \emph{Psychological Methods},
#' \bold{10}, 40-64.
#'
#' Cheung, M. W.-L., & Chan, W. (2009). A two-stage approach to synthesizing
#' covariance matrices in meta-analytic structural equation modeling.
#' \emph{Structural Equation Modeling}, \bold{16}, 28-53.
#' @source Inter-University Consortium for Political and Social Research.
#' (1989). \emph{International Social Survey Program: Work orientation}. Ann
#' Arbor, MI: Author.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(issp89)
#'
#' #### Analysis of correlation structure in Cheung and Chan (2005)
#' #### Fixed-effects model: Stage 1 analysis
#' cor1 <- tssem1(issp89$data, issp89$n, method="FEM", cor.analysis=TRUE)
#' summary(cor1)
#'
#' ## Three-factor CFA model in lavaan syntax
#' ## f1: Job Prospects (JP1-JP3); f2: Job Nature (JN1-JN4); f3: Time Demand (TD1-TD2)
#' model <- "## Factor loadings
#'           f1 =~ f1JP1*JP1 + f1JP2*JP2 + f1JP3*JP3
#'           f2 =~ f2JN1*JN1 + f2JN2*JN2 + f2JN3*JN3 + f2JN4*JN4
#'           f3 =~ f3TD1*TD1 + f3TD2*TD2
#'           ## Factor correlations
#'           f1 ~~ corf2f1*f2
#'           f1 ~~ corf3f1*f3
#'           f2 ~~ corf3f2*f3"
#'
#' RAM1 <- lavaan2RAM(model, obs.variables=c("JP1","JP2","JP3",
#'                                           "JN1","JN2","JN3","JN4",
#'                                           "TD1","TD2"))
#' RAM1
#'
#' ## ## Equivalent RAM specification using create.mxMatrix()
#' ## Phi <- create.mxMatrix( c("0.3*corf2f1","0.3*corf3f1","0.3*corf3f2"),
#' ##                         type="Stand", as.mxMatrix=FALSE )
#' ## Psi <- create.mxMatrix( paste("0.2*e", 1:9, sep=""), type="Diag",
#' ##                         as.mxMatrix=FALSE )
#' ## S1 <- bdiagMat(list(Psi, Phi))
#' ## S1 <- as.mxMatrix(S1)
#' ##
#' ## Lambda <- create.mxMatrix( c(".3*f1x1",".3*f1x2",".3*f1x3",rep(0,9),
#' ##                              ".3*f2x4",".3*f2x5",".3*f2x6",".3*f2x7",
#' ##                              rep(0,9),".3*f3x8",".3*f3x9"), type="Full",
#' ##                              ncol=3, nrow=9, as.mxMatrix=FALSE )
#' ## Zero1 <- matrix(0, nrow=9, ncol=9)
#' ## Zero2 <- matrix(0, nrow=3, ncol=12)
#' ## A1 <- rbind( cbind(Zero1, Lambda), Zero2 )
#' ## A1 <- as.mxMatrix(A1)
#' ##
#' ## F1 <- create.Fmatrix(c(rep(1,9), rep(0,3)))
#'
#' #### Fixed-effects model: Stage 2 analysis
#' cor2 <- tssem2(cor1, RAM=RAM1, intervals.type="LB")
#' summary(cor2)
#'
#' ## Display the model with the parameter estimates
#' plot(cor2, nDigits=1)
#'
#' #### Analysis of covariance structure in Cheung and Chan (2009)
#' #### Fixed-effects model: Stage 1 analysis
#' cov1 <- tssem1(issp89$data, issp89$n, method="FEM", cor.analysis=FALSE)
#' summary(cov1)
#'
#' #### Fixed-effects model: Stage 2 analysis
#' cov2 <- tssem2(cov1, RAM=RAM1)
#' summary(cov2)
#'
#' ## Display the model with the parameter estimates
#' plot(cov2, nDigits=1)
#' }
#'
NULL





#' Dataset from Jaramillo, Mulki and Marshall (2005)
#'
#' A dataset of the relationship between organizational commitment (OC) and
#' salesperson job performance (JP) from Jaramillo, Mulki & Marshall (2005).
#'
#'
#' @name Jaramillo05
#' @docType data
#' @format A data frame with 61 observations on the following 10 variables.
#' \describe{ \item{\code{Author}}{a character vector of study}
#' \item{\code{Sample_size}}{sample size of the study}
#' \item{\code{Sales}}{sample type; either "mixed", "nonsales" or "sales"}
#' \item{\code{Country}}{a character vector of country of study}
#' \item{\code{IDV}}{Hofstede's (1997) individualism index}
#' \item{\code{OC_scale}}{scale of OC; either "Porter or Mowday", "Meyer" or
#' "other"} \item{\code{OC_alpha}}{Coefficient alpha of organizational
#' commitment} \item{\code{JP_alpha}}{Coefficient alpha of job performance}
#' \item{\code{r}}{correlation between organizational commitment and job
#' performance} \item{\code{r_v}}{sampling variance of r}
#' \item{Citations}{Citations from Google Scholar as of 27 August 2024} }
#' @source Jaramillo, F., Mulki, J. P., & Marshall, G. W. (2005). A
#' meta-analysis of the relationship between organizational commitment and
#' salesperson job performance: 25 years of research. \emph{Journal of Business
#' Research}, \bold{58(6)}, 705-714. doi:10.1016/j.jbusres.2003.10.004
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' ## Research question 4.4.1
#' summary(meta(r, r_v, data=Jaramillo05))
#'
#' ## Research question 4.4.2
#' ## Select cases with either "sales" or "nonsales"
#' Sales.df <- subset(Jaramillo05, Sales %in% c("sales", "nonsales"))
#'
#' ## Create a predictor with 1 and 0 when they are "sales" or "nonsales", respectively
#' predictor <- ifelse(Jaramillo05$Sales=="sales", yes=1, no=0)
#'
#' ## Mixed-effects meta-analysis
#' summary( meta(y = r, v = r_v, x = predictor, data = Jaramillo05) )
#'
#' ## Research question 4.4.3
#' summary(meta(r, r_v, x=IDV, data=Jaramillo05))
#' }
#'
NULL





#' Multivariate effect sizes reported by Kalaian and Raudenbush (1996)
#'
#' This data set includes 47 multivariate effect sizes reported by Kalaian and
#' Raudenbush (1996, Table 1).
#'
#' A list of data with the following structure: \describe{ \item{Study}{Study
#' name} \item{Year}{Year of publication} \item{n_e}{Sample size of the
#' experimental group} \item{n_c}{Sample size of the control group}
#' \item{dSAT_V}{Standardized mean difference of the Scholastic Aptitude Test
#' (SAT) on verbal} \item{dSAT_M}{Standardized mean difference of SAT on math}
#' \item{var_V}{Sampling variance of \code{dSAT_V}} \item{cov_VM}{Sampling
#' covariance of \code{dSAT_V} and \code{dSAT_M} with a common correlation of
#' 0.66} \item{var_M}{Sampling variance of \code{dSAT_M}} \item{Hr}{Hours of
#' training} \item{ETS}{Educational Testing Service} \item{Study_type}{Either
#' \code{Randomized}, \code{Matched} or \code{Nonequivalent comparison}}
#' \item{Home_work}{Home work} }
#'
#' @name Kalaian96
#' @docType data
#' @source Kalaian, H. A., & Raudenbush, S. W. (1996). A multivariate mixed
#' linear model for meta-analysis. \emph{Psychological Methods}, \emph{1}(3),
#' 227-235. https://doi.org/10.1037/1082-989X.1.3.227
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Kalaian96)
#' }
#'
NULL





#' Eight studies from Mak et al. (2009)
#'
#' Eight studies from Mak et al. (2009) were reported by Cheung et al. (2012).
#'
#'
#' @name Mak09
#' @docType data
#' @format A data frame with 8 observations on the following 10 variables.
#' \describe{ \item{\code{Study}}{a character vector of study}
#' \item{\code{type}}{a character vector} \item{\code{AF.BP}}{a numeric
#' vector} \item{\code{Tot.BP}}{a numeric vector} \item{\code{AF.non.BP}}{a
#' numeric vector} \item{\code{Tot.non.BP}}{a numeric vector}
#' \item{\code{yi}}{a numeric vector} \item{\code{vi}}{a numeric vector}
#' \item{\code{age.mean}}{a numeric vector} \item{\code{study.duration}}{a
#' numeric vector} }
#' @references Cheung, M. W.-L., Ho, R. C. M., Lim, Y., & Mak, A. (2012).
#' Conducting a meta-analysis: Basics and good practices. \emph{International
#' Journal of Rheumatic Diseases}, \bold{15(2)}, 129-135. doi:
#' 10.1111/j.1756-185X.2012.01712.x
#' @source Mak, A., Cheung, M. W.-L., Ho, R. C. M., Cheak, A. A. C., & Lau, C.
#' S. (2009). Bisphosphonate and atrial fibrillation: Bayesian meta-analyses of
#' randomized controlled trials and observational studies. \emph{BMC
#' Musculoskeletal Disorders}, \bold{10(113)}. doi:10.1186/1471-2474-10-113
#' Available at
#' \url{https://link.springer.com/article/10.1186/1471-2474-10-113}.
#' @keywords datasets
#' @examples
#'
#' ## Random-effects meta-analysis
#' ( meta1 <- summary(meta(y=yi, v=vi, data=Mak09, I2=c("I2q", "I2hm"))) )
#'
#' ## Convert the estimates back into odds ratio 
#' OR <- with(coef(meta1), exp(c(Estimate[1], lbound[1], ubound[1])))
#' names(OR) <- c("Estimate in OR", "lbound in OR", "ubound in OR")
#' OR
#'
#' ## Mixed-effects meta-analysis with mean age as a predictor
#' summary( meta(y=yi, v=vi, x=age.mean, data=Mak09) )
#'
NULL





#' Correlation Matrices from Mathieu et al. (2015)
#'
#' The data set includes a list of correlation matrices of panel studies
#' between cohesion (C) and performance (P) in Mathieu et al. (2015, Table 1).
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of studies of correlation matrices. The variables are \emph{C1}, \emph{P1},
#' \emph{C2}, and \emph{P2}.} \item{n}{A vector of sample sizes.}
#' \item{Year}{Year of publication.} \item{Sample}{Sample characteristics.}
#' \item{Student}{Whether the samples are student or non-student based on
#' \code{Sample}.} }
#'
#' @name Mathieu15
#' @docType data
#' @source Mathieu, J. E., Kukenberger, M. R., D'Innocenzo, L., & Reilly, G.
#' (2015). Modeling reciprocal team cohesion-performance relationships, as
#' impacted by shared leadership and members' competence. \emph{Journal of
#' Applied Psychology}, \bold{100}(3), 713-734.
#' https://doi.org/10.1037/a0038898
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' # TSSEM
#' ## Model 1: no constraint
#' ## Stage 1 analysis
#' tssem1.fit <- tssem1(Mathieu15$data, Mathieu15$n)
#' summary(tssem1.fit)
#'
#' ## Proposed model in lavaan syntax
#' model1 <- 'C2 ~ c2c*C1 + p2c*P1
#'            P2 ~ c2p*C1 + p2p*P1
#'            C1 ~~ c1withp1*P1
#'            C1 ~~ 1*C1
#'            P1 ~~ 1*P1
#'            C2 ~~ c2withp2*P2'
#'
#' ## Convert the lavaan model to RAM specification
#' RAM1 <- lavaan2RAM(model1, obs.variables=c("C1", "P1", "C2", "P2"))
#' RAM1
#'
#' ## Stage 2 analysis
#' tssem1b.fit <- tssem2(tssem1.fit, RAM=RAM1)
#' summary(tssem1b.fit)
#'
#' plot(tssem1b.fit, col="yellow", edge.label.position=0.58)
#'
#' ## Model 2: Equality constraints on the path coefficient
#' ## Proposed model with equal effects time 1 to time 2
#' model2 <- 'C2 ~ same*C1 + diff*P1
#'            P2 ~ diff*C1 + same*P1
#'            C1 ~~ c1withp1*P1
#'            C1 ~~ 1*C1
#'            P1 ~~ 1*P1
#'            C2 ~~ c2withp2*P2'
#'
#' ## Convert the lavaan model to RAM specification
#' RAM2 <- lavaan2RAM(model2, obs.variables=c("C1", "P1", "C2", "P2"))
#' RAM2
#'
#' ## Stage 2 analysis
#' tssem2b.fit <- tssem2(tssem1.fit, RAM=RAM2)
#' summary(tssem2b.fit)
#'
#' ## Compare the models with and without the constraints. 
#' anova(tssem1b.fit, tssem2b.fit)
#'
#' ## Plot the model
#' plot(tssem2b.fit, col="yellow", edge.label.position=0.60)
#'
#'
#' ## OSMASEM
#' my.df <- Cor2DataFrame(Mathieu15)
#'
#' head(my.df$data)
#'
#' ## Model without any moderator
#' osmasem.fit1 <- osmasem(model.name="No moderator", RAM=RAM1, data=my.df)
#' summary(osmasem.fit1)
#'
#' ## Extract the heterogeneity variance-covariance matrix
#' diag(VarCorr(osmasem.fit1))
#'
#' plot(osmasem.fit1, col="yellow", edge.label.position=0.6)
#'
#' ## Model with student sample as a moderator on the regression coefficients
#' A1 <- create.modMatrix(RAM1, output="A", "Student")
#' A1
#'
#' ## Model with a moderator    
#' osmasem.fit2 <- osmasem(model.name="Student sample as a moderator", RAM=RAM1, 
#'                         Ax=A1, data=my.df)
#' summary(osmasem.fit2)
#'
#' ## Compare the models with and without the moderator
#' anova(osmasem.fit2, osmasem.fit1)
#'
#' ## Get the R2 of the moderator
#' osmasemR2(osmasem.fit2, osmasem.fit1)
#' }
#'
NULL





#' Meta-Analysis using Structural Equation Modeling
#'
#' A collection of functions for conducting meta-analysis using a structural
#' equation modeling (SEM) approach via the 'OpenMx' and 'lavaan' packages. It
#' also implements various procedures to perform meta-analytic structural
#' equation modeling on the correlation and covariance matrices.
#'
#' \tabular{ll}{ Package: \tab metaSEM\cr Type: \tab Package\cr Version: \tab
#' 1.5.5\cr Date: \tab 2026-05-02\cr License: \tab GPL (>=2)\cr LazyLoad: \tab
#' yes\cr }
#'
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#'
#' Maintainer: Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @references Cheung, M. W.-L. (2008). A model for integrating fixed-,
#' random-, and mixed-effects meta-analyses into structural equation modeling.
#' \emph{Psychological Methods}, \bold{13} (3), 182-202.
#' https://doi.org/10.1037/a0013163
#'
#' Cheung, M. W.-L. (2009). Constructing approximate confidence intervals for
#' parameters with structural equation models. \emph{Structural Equation
#' Modeling}, \bold{16} (2), 267-294. https://doi.org/10.1080/10705510902751291
#'
#' Cheung, M. W.-L. (2010). Fixed-effects meta-analyses as multiple-group
#' structural equation models. \emph{Structural Equation Modeling}, \bold{17}
#' (3), 481-509. https://doi.org/10.1080/10705511.2010.489367
#'
#' Cheung, M. W.-L. (2013). Implementing restricted maximum likelihood
#' estimation in structural equation models. \emph{Structural Equation
#' Modeling}, \bold{20} (1), 157-167.
#' https://doi.org/10.1080/10705511.2013.742404
#'
#' Cheung, M. W.-L. (2013). Multivariate meta-analysis as structural equation
#' models. \emph{Structural Equation Modeling}, \bold{20} (3), 429-454.
#' https://doi.org/10.1080/10705511.2013.797827
#'
#' Cheung, M. W.-L. (2014). Modeling dependent effect sizes with three-level
#' meta-analyses: A structural equation modeling approach. \emph{Psychological
#' Methods}, \bold{19} (2), 211-229. https://doi.org/10.1037/a0032968
#'
#' Cheung, M. W.-L. (2014). Fixed- and random-effects meta-analytic structural
#' equation modeling: Examples and analyses in R. \emph{Behavior Research
#' Methods}, \bold{46} (1), 29-40. https://doi.org/10.3758/s13428-013-0361-y
#'
#' Cheung, M. W.-L. (2015). metaSEM: An R package for meta-analysis using
#' structural equation modeling. \emph{Frontiers in Psychology}, \bold{5}
#' (1521). https://doi.org/10.3389/fpsyg.2014.01521
#'
#' Cheung, M. W.-L. (2015). \emph{Meta-Analysis: A Structural Equation Modeling
#' Approach}. Chichester, West Sussex: John Wiley & Sons, Inc.
#'
#' Cheung, M. W.-L. (2018). Issues in solving the problem of effect size
#' heterogeneity in meta-analytic structural equation modeling: A commentary
#' and simulation study on Yu, Downes, Carter, and O'Boyle (2016).
#' \emph{Journal of Applied Psychology}, \bold{103} (7), 787-803.
#' https://doi.org/10.1037/apl0000284
#'
#' Cheung, M. W.-L. (2018). Computing multivariate effect sizes and their
#' sampling covariance matrices with structural equation modeling: Theory,
#' examples, and computer simulations. \emph{Frontiers in Psychology}, \bold{9}
#' (1387). https://doi.org/10.3389/fpsyg.2018.01387
#'
#' Cheung, M. W.-L. (2019). Some reflections on combining meta-analysis and
#' structural equation modeling. \emph{Research Synthesis Methods}, \bold{10}
#' (1), 15-22. https://doi.org/10.1002/jrsm.1321
#'
#' Cheung, M. W.-L. (2021). Meta-analytic structural equation modeling. In
#' \emph{Oxford Research Encyclopedia of Business and Management}. Oxford
#' University Press. https://doi.org/10.1093/acrefore/9780190224851.013.225
#'
#' Cheung, M. W.-L., & Chan, W. (2004). Testing dependent correlation
#' coefficients via structural equation modeling. \emph{Organizational Research
#' Methods}, \bold{7} (2), 206-223. https://doi.org/10.1177/1094428104264024
#'
#' Cheung, M. W.-L., & Chan, W. (2005). Meta-analytic structural equation
#' modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10} (1),
#' 40-64. https://doi.org/10.1037/1082-989X.10.1.40
#'
#' Cheung, M. W.-L., & Chan, W. (2009). A two-stage approach to synthesizing
#' covariance matrices in meta-analytic structural equation modeling.
#' \emph{Structural Equation Modeling}, \bold{16} (1), 28-53.
#' https://doi.org/10.1080/10705510802561295
#'
#' Cheung, M. W.-L., & Cheung, S.-F. (2016). Random-effects models for
#' meta-analytic structural equation modeling: Review, issues, and
#' illustrations. \emph{Research Synthesis Methods}, \bold{7} (2), 140-155.
#' https://doi.org/10.1002/jrsm.1166
#'
#' Jak, S., & Cheung, M. W.-L. (2018). Testing moderator hypotheses in
#' meta-analytic structural equation modeling using subgroup analysis.
#' \emph{Behavior Research Methods}, \bold{50} (4), 1359-1373.
#' https://doi.org/10.3758/s13428-018-1046-3
#'
#' Jak, S., & Cheung, M. W.-L. (2020). Meta-analytic structural equation
#' modeling with moderating effects on SEM parameters. \emph{Psychological
#' Methods}, \bold{25} (4), 430-455. https://doi.org/10.1037/met0000245
#'
#' Jak, S., & Cheung, M. W.-L. (2024). A cautionary note on using univariate
#' methods for meta-analytic structural equation modeling. \emph{Advances in
#' Methods and Practices in Psychological Science}, \bold{7}(4), 1-24.
#' https://doi.org/10.1177/25152459241274249
#' @keywords internal
"_PACKAGE"





#' Dataset on the Environmental Tobacco Smoke (ETS) on children's health
#'
#' This dataset includes 59 studies reported by Nam, Mengersen, and Garthwaite
#' (2003) on the potential health effects among children exposed to
#' environmental tobacco smoke (ETS), or passive smoking. The effect sizes are
#' the log odds ratios of asthma and lower respiratory disease (LRD).
#'
#' A list of data with the following structure: \describe{ \item{ID}{Study
#' identification number.} \item{Size}{Total number of valid subjects in the
#' study.} \item{Age}{Mean age of participants.} \item{Year}{Year of
#' publication.} \item{Country}{Country code.} \item{Smoke}{Source of ETS.}
#' \item{Adj}{Whether the reported odds ratio is adjusted for covariates.}
#' \item{Asthma_logOR}{Log odds ratio of asthma.} \item{LRD_logOR}{Log odds
#' ratio of lower respiratory disease.} \item{Asthma_v}{Sampling variance of
#' Asthma_logOR.} \item{AsthmaLRD_cov_05}{Sampling covariance between
#' Asthma_logOR and LRD_logOR by assuming a correlation of 0.5}
#' \item{LRD_v}{Sampling variance of LRD_logOR.} }
#'
#' @name Nam03
#' @docType data
#' @source Nam, I.-S., Mengersen, K., & Garthwaite, P. (2003). Multivariate
#' meta-analysis. \emph{Statistics in Medicine}, \bold{22}(14), 2309-2333.
#' https://doi.org/10.1002/sim.1410
#' @keywords datasets
#' @examples
#'
#' data(Nam03)
#'
NULL





#' Correlation Matrices from Nohe et al. (2015)
#'
#' The data sets include two lists of correlation matrices of panel studies
#' between work-family conflict and strain reported in Table A1
#' (\code{Nohe15A1}) and Table A2 (\code{Nohe15A2}) by Nohe et al. (2015).
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of studies of correlation matrices. The variables are \emph{W1}, \emph{S1},
#' \emph{W2}, and \emph{S2} in \code{Nohe15A1} and \emph{F1}, \emph{S1},
#' \emph{F2}, and \emph{S2} in \code{Nohe15A2}} \item{n}{A vector of sample
#' sizes} \item{RelXX}{The reliabilities of \emph{W1}, \emph{S1}, \emph{W2} and
#' \emph{S2} in \code{Nohe15A1} and the reliabilities of \emph{F1} \emph{S1},
#' \emph{F2} , and \emph{S2} in \code{Nohe15A2}} \item{FemalePer}{Percentage of
#' female participants} \item{Publication}{Whether the studies were published
#' (\emph{P}) or unpublished (\emph{U})} \item{Lag}{Time lag between the coded
#' measurement waves in months} }
#'
#' @name Nohe15
#' @aliases Nohe15A1 Nohe15A2
#' @docType data
#' @source Nohe, C., Meier, L. L., Sonntag, K., & Michel, A. (2015). The
#' chicken or the egg? A meta-analysis of panel studies of the relationship
#' between work-family conflict and strain. \emph{Journal of Applied
#' Psychology}, \bold{100}(2), 522-536.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' #### TSSEM
#'
#' ## Set seed for replicability    
#' set.seed(23891)
#'
#' ## Table A1
#' randA1a <- tssem1(Nohe15A1$data, Nohe15A1$n, method="REM", RE.type="Diag")
#' summary(randA1a)
#'
#' model1 <- 'W2 ~ w2w*W1 + s2w*S1
#'            S2 ~ w2s*W1 + s2s*S1
#'            W1 ~~ w1WITHs1*S1
#'            W2 ~~ w2WITHs2*S2
#'            W1 ~~ 1*W1
#'            S1 ~~ 1*S1
#'            W2 ~~ Errw2*W2
#'            S2 ~~ Errs2*S2'
#'
#' ## Display the model
#' plot(model1, layout="spring")    
#'
#' RAM1 <- lavaan2RAM(model1, obs.variables=c("W1", "S1", "W2", "S2"))
#' RAM1
#'
#' randA1b <- tssem2(randA1a, RAM=RAM1)
#' summary(randA1b)
#'
#' ## Display the model with the parameter estimates
#' plot(randA1b, layout="spring")    
#'
#' ## Table A2
#' randA2a <- tssem1(Nohe15A2$data, Nohe15A2$n, method="REM", RE.type="Diag")
#' ## Rerun to remove error code
#' randA2a <- rerun(randA2a)
#' summary(randA2a)
#'
#' model2 <- 'F2 ~ f2f*F1 + s2F*S1
#'            S2 ~ f2s*F1 + s2s*S1
#'            F1 ~~ f1WITHs1*S1
#'            F2 ~~ f2WITHs2*S2
#'            F1 ~~ 1*F1
#'            S1 ~~ 1*S1
#'            F2 ~~ Errf2*F2
#'            S2 ~~ Errs2*S2'
#'
#' ## Display the model
#' plot(model2, layout="spring")
#'
#' RAM2 <- lavaan2RAM(model2, obs.variables=c("F1", "S1", "F2", "S2"))
#' RAM2
#'
#' randA2b <- tssem2(randA2a, RAM=RAM2)
#' summary(randA2b)
#'
#' ## Display the model with the parameter estimates
#' plot(randA2b, layout="spring")  
#'
#' ## Estimate the heterogeneity of the parameter estimates
#' tssemParaVar(randA1a, randA2b)    
#'
#' ## Parametric bootstrap based on Yu et al. (2016)
#' ## I assume that you know what you are doing!
#'
#' ## Set seed for reproducibility
#' set.seed(39128482)
#'
#' ## Average the correlation coefficients with the univariate-r approach
#' uni1 <- uniR1(Nohe15A1$data, Nohe15A1$n)
#' uni1
#'
#' ## Generate random correlation matrices
#' boot.cor <- bootuniR1(uni1, Rep=50)
#'
#' ## Display the quality of the generated correlation matrices
#' summary(boot.cor)
#'
#' ## Proposed saturated model
#' model1 <- 'W2 + S2 ~ W1 + S1'
#'
#' ## Use the harmonic mean of the sample sizes as n in SEM
#' n <- uni1$n.harmonic    
#'
#' boot.fit1 <- bootuniR2(model=model1, data=boot.cor, n=n)
#' summary(boot.fit1)
#'
#' ## Proposed model with equal regression coefficients
#' model2 <- 'W2 ~ Same*W1 + Cross*S1
#'            S2 ~ Cross*W1 + Same*S1'
#'
#' boot.fit2 <- bootuniR2(model=model2, data=boot.cor, n=n)
#' summary(boot.fit2)
#'
#' #### OSMASEM    
#'
#' ## Calculate the sampling variance-covariance matrix of the correlation matrices.    
#' my.df <- Cor2DataFrame(Nohe15A1)
#'
#' ## Standardize the moderator "Lag"
#' my.df$data$Lag <- scale(my.df$data$Lag)
#'
#' head(my.df$data)
#'
#' ## Proposed model
#' model1 <- 'W2 ~ w2w*W1 + s2w*S1
#'            S2 ~ w2s*W1 + s2s*S1
#'            W1 ~~ w1WITHs1*S1
#'            W2 ~~ w2WITHs2*S2
#'            W1 ~~ 1*W1
#'            S1 ~~ 1*S1
#'            W2 ~~ Errw2*W2
#'            S2 ~~ Errs2*S2'
#' plot(model1)     
#'
#' ## Convert it into RAM specification    
#' RAM1 <- lavaan2RAM(model1, obs.variables=c("W1", "S1", "W2", "S2"))
#' RAM1
#'
#' ## Create vechs of the model implied correlation matrix
#' ## with implicit diagonal constraints
#' ## M0 <- create.vechsR(A0=RAM1$A, S0=RAM1$S)
#'
#' ## Create heterogeneity variances
#' ## RE.type= either "Diag" or "Symm"
#' ##
#' ## Transform= either "expLog" or "sqSD" for better estimation on variances
#' ## T0 <- create.Tau2(RAM=RAM1, RE.type="Diag")
#' ##
#' ## Fit the model    
#' ## fit0 <- osmasem(model.name="No moderator", Mmatrix=M0, Tmatrix=T0, data=my.df)
#'
#' ## Fit the model
#' fit0 <- osmasem(model.name="No moderator", RAM=RAM1, data=my.df)
#' summary(fit0)
#'
#' ## Get the SRMR
#' osmasemSRMR(fit0)
#'
#' ## Get the transformed variance component of the random effects    
#' VarCorr(fit0)
#'
#' ## "lag" as a moderator on A matrix
#' A1 <- matrix(c(0,0,0,0,
#'                0,0,0,0,
#'                "0*data.Lag","0*data.Lag",0,0,
#'                "0*data.Lag","0*data.Lag",0,0),
#'              nrow=4, ncol=4, byrow=TRUE)
#'
#' ## M1 <- create.vechsR(A0=RAM1$A, S0=RAM1$S, Ax=A1)
#' ##
#' ## Fit the nodel
#' ## fit1 <- osmasem(model.name="Lag as a moderator for Amatrix", Mmatrix=M1,
#' ##                 Tmatrix=T0, data= my.df)
#'
#' fit1 <- osmasem(model.name="Lag as a moderator for Amatrix",
#'                 RAM=RAM1, Ax=A1, data= my.df)
#' summary(fit1)
#' VarCorr(fit1)
#'
#' ## Compare the models with and without the moderator "lag"
#' anova(fit1, fit0)
#'
#' ## Calculate the R2    
#' osmasemR2(fit0, fit1)
#' }
#'
NULL





#' Studies on the Hospital Anxiety and Depression Scale Reported by Norton et
#' al. (2013)
#'
#' The data set includes 28 studies on 14 items measuring the Hospital Anxiety
#' and Depression Scale (HADS) Reported by Norton et al. (2013).
#'
#' The variables are: \describe{ \item{data}{A list of 28 studies of
#' correlation matrices. The variables are 14 items (x1 to x14) measuring
#' HADS.} \item{n}{A vector of sample sizes} \item{population}{A vector of the
#' population of the data} \item{group}{A vector of classification into
#' \emph{patients} vs. \emph{non-patients} based on population} }
#'
#' @name Norton13
#' @docType data
#' @references Jak, S., & Cheung, M. W.-L. (2018). Addressing heterogeneity in
#' meta-analytic structural equation modeling using subgroup analysis.
#' \emph{Behavior Research Methods}, \bold{50}, 1359-1373.
#' @source Norton, S., Cosco, T., Doyle, F., Done, J., & Sacker, A. (2013). The
#' Hospital Anxiety and Depression Scale: A meta confirmatory factor analysis.
#' \emph{Journal of Psychosomatic Research}, \emph{74}(1), 74-81.
#' @keywords datasets
#' @examples
#'
#' data(Norton13)
#'
NULL





#' Plot methods for various objects
#'
#' It plots the models from either the lavaan model or \code{meta}, \code{wls},
#' and \code{osmasem} objects.
#'
#'
#' @name plot
#' @aliases plot.meta plot.character plot.wls plot.osmasem plot.osmasem2 plot.mxsem
#' @param x An object returned from either a lavaan model class
#' \code{character}, \code{osmasem}, \code{osmasem3L}, \code{wls} or
#' \code{meta}
#' @param effect.sizes Numeric values indicating which effect sizes to be
#' plotted. At least two effect sizes are required. To plot the effect sizes of
#' \eqn{y_1}{y1} and \eqn{y_2}{y2}, one may use \code{effect.sizes=c(1,2)}. If
#' it is missing, all effect sizes will be plotted in a pairwise way.
#' @param add.margin Value for additional margins on the left and bottom
#' margins.
#' @param interval Interval for the confidence ellipses.
#' @param main Main title of each plot. If there are multiple plots, a vector
#' of character titles may be used.
#' @param axis.labels Labels for the effect sizes.
#' @param study.col The color for individual studies. See \code{col} in
#' \code{\link[graphics]{par}}.
#' @param study.pch Plotting character of individual studies. See \code{pch} in
#' \code{\link[graphics]{points}}.
#' @param study.min.cex The minimum value of cex for individual studies. See
#' \code{cex} in \code{\link[graphics]{par}}.
#' @param study.weight.plot Logical. If \code{TRUE}, the plotting size of
#' individual studies (cex) will be proportional to one over the square root of
#' the determinant of the sampling covariance matrix of the study.
#' @param study.ellipse.plot Logical. If \code{TRUE}, the confidence ellipses
#' of individual studies are plotted.
#' @param study.ellipse.col The color of the confidence ellipses of individual
#' studies. See \code{col} in \code{\link[graphics]{par}}.
#' @param study.ellipse.lty The line type of the confidence ellipse of
#' individual studies. See \code{lty} in \code{\link[graphics]{par}}.
#' @param study.ellipse.lwd The line width of the confidence ellipse of
#' individual studies. See \code{lwd} in \code{\link[graphics]{par}}.
#' @param estimate.col The color of the estimated effect size. See \code{col}
#' in \code{\link[graphics]{par}}.
#' @param estimate.pch Plotting character of the estimated effect sizes. See
#' \code{pch} in \code{\link[graphics]{points}}.
#' @param estimate.cex The amount of plotting of the estimated effect sizes.
#' See \code{cex} in \code{\link[graphics]{par}}.
#' @param estimate.ellipse.plot Logical. If \code{TRUE}, the confidence ellipse
#' of the estimated effect sizes will be plotted.
#' @param estimate.ellipse.col The color of the confidence ellipse of the
#' estimated effect sizes. See \code{col} in \code{\link[graphics]{par}}.
#' @param estimate.ellipse.lty The line type of the confidence ellipse of the
#' estimated effect sizes. See \code{lty} in \code{\link[graphics]{par}}.
#' @param estimate.ellipse.lwd The line width of the confidence ellipse of the
#' estimated effect sizes. See \code{lwd} in \code{\link[graphics]{par}}.
#' @param randeff.ellipse.plot Logical. If \code{TRUE}, the confidence ellipses
#' of the random effects will be plotted.
#' @param randeff.ellipse.col Color of the confidence ellipses of the random
#' effects. See \code{col} in \code{\link[graphics]{par}}.
#' @param randeff.ellipse.lty The line type of the confidence ellipses of the
#' random effects. See \code{lty} in \code{\link[graphics]{par}}.
#' @param randeff.ellipse.lwd The line width of the confidence ellipses of the
#' random effects. See \code{lwd} in \code{\link[graphics]{par}}.
#' @param univariate.plot Logical. If \code{TRUE}, the estimated univariate
#' effect sizes will be plotted.
#' @param univariate.lines.col The color of the estimated univariate effect
#' sizes. See \code{col} in \code{\link[graphics]{par}}.
#' @param univariate.lines.lty The line type of the estimated univariate effect
#' sizes. See \code{lty} in \code{\link[graphics]{par}}.
#' @param univariate.lines.lwd The line width of the estimated univariate
#' effect sizes. See \code{lwd} in \code{\link[graphics]{par}}.
#' @param univariate.polygon.width The width of the polygon of the estimated
#' univariate effect sizes.
#' @param univariate.polygon.col The color of the polygon of the estimated
#' univariate effect sizes.
#' @param univariate.arrows.col The color of the arrows of the estimated
#' univariate effect sizes.
#' @param univariate.arrows.lwd The line width of the arrows of the estimated
#' univariate effect sizes.
#' @param diag.panel Logical. If \code{TRUE}, diagonal panels will be created.
#' They can then be used for forrest plots for univariate meta-analysis.
#' @param xlim NULL or a numeric vector of length 2; if it is NULL, it provides
#' defaults estimated from the data.
#' @param ylim NULL or a numeric vector of length 2; if it is NULL, it provides
#' defaults estimated from the data.
#' @param fixed.x Argument passed to \code{\link[semPlot]{semPlotModel}}.
#' @param manNames Argument passed to \code{\link[semPlot]{semPaths}}
#' @param latNames Argument passed to \code{\link[semPlot]{semPaths}}
#' @param labels Argument passed to \code{\link[semPlot]{semPaths}}
#' @param what Argument passed to \code{\link[semPlot]{semPaths}}
#' @param nCharNodes Argument passed to \code{\link[semPlot]{semPaths}}
#' @param nCharEdges Argument passed to \code{\link[semPlot]{semPaths}}
#' @param layout Argument passed to \code{\link[semPlot]{semPaths}}
#' @param color Argument passed to \code{\link[semPlot]{semPaths}}
#' @param sizeMan Argument passed to \code{\link[semPlot]{semPaths}}
#' @param sizeLat Argument passed to \code{\link[semPlot]{semPaths}}
#' @param edge.label.cex Argument passed to \code{\link[semPlot]{semPaths}}
#' @param weighted Argument passed to \code{\link[semPlot]{semPaths}}
#' @param \dots Further arguments passed to the methods.
#' @note The estimated effect sizes and random effects are based on the labels
#' Intercept1, Intercept2, ... and Tau2_1_1, Tau2_2_1, Tau2_2_2, etc. At least
#' two effect sizes are required for this function.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{Berkey98}}, \code{\link[metaSEM]{wvs94a}}
#' \code{\link[metaSEM]{meta2semPlot}} \code{\link[semPlot]{semPaths}}
#' @references Cheung, M. W.-L. (2013). Multivariate meta-analysis as
#' structural equation models. \emph{Structural Equation Modeling}, \bold{20},
#' 429-454.
#' @keywords methods
#' @examples
#'
#' \donttest{
#' ## lavaan model
#' model <- "y ~ m + x
#'           m ~ x"
#' plot(model)
#' }
#'
NULL





#' Print Methods for various Objects
#'
#' Print methods for the \code{tssem1FEM}, \code{tssem1FEM.cluster},
#' \code{tssem1REM}, \code{wls}, \code{meta}, \code{meta3LFIML}, \code{reml},
#' \code{uniR1} and \code{impliedR} objects.
#'
#'
#' @name print
#' @aliases print.tssem1FEM print.tssem1FEM.cluster print.tssem1REM print.wls print.meta print.meta3LFIML print.reml print.uniR1 print.impliedR
#' @param x An object returned from either class \code{tssem1FEM}, class
#' \code{tssem1FEM.cluster}, class \code{tssem1REM}, class \code{wls}, class
#' \code{meta}, class \code{meta3LFIML}, class \code{reml}, class \code{uniR1}
#' or class \code{impliedR}
#' @param \dots Further arguments to be passed to \code{summary.default} or
#' unused.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
#' \code{\link[metaSEM]{meta}}, \code{\link[metaSEM]{reml}}
#' @keywords methods
NULL





#' Read External Correlation/Covariance Matrices
#'
#' It reads full/lower triangle/stacked vectors of correlation/covariance data
#' into a list of correlation/covariance matrices.
#'
#'
#' @name readData
#' @aliases readFullMat readStackVec readLowTriMat
#' @param file File name of the data.
#' @param no.var The number of variables in the data.
#' @param \dots Further arguments to be passed to \code{\link[base]{scan}} for
#' \code{readLowTriMat} and to \code{\link[utils]{read.table}} for
#' \code{readFullMat} and \code{readStackVec}.
#' @return A list of correlation/covariance matrices.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @keywords utilities
#' @examples
#'
#' \donttest{
#' ## Write two full correlation matrices into a file named "fullmat.dat".
#' ## x2 is missing in the second matrix.
#' ## The content of "fullmat.dat" is
#' # 1.0 0.3 0.4
#' # 0.3 1.0 0.5
#' # 0.4 0.5 1.0
#' # 1.0 NA 0.4
#' # NA NA NA
#' # 0.4 NA 1.0
#'
#' ## cat("1.0 0.3 0.4\n0.3 1.0 0.5\n0.4 0.5 1.0
#' ## 1.0 NA 0.4\nNA NA NA\n0.4 NA 1.0",
#' ## file="fullmat.dat", sep="")
#'
#' ## Read the correlation matrices from a file
#' ## my.full <- readFullMat("fullmat.dat")
#'
#' ## Read the correlation matrices from a string
#' x <-
#' "1.0 0.3 0.4
#' 0.3 1.0 0.5
#' 0.4 0.5 1.0
#' 1.0 NA 0.4
#' NA NA NA
#' 0.4 NA 1.0"
#'
#' my.full <- readFullMat(textConnection(x))
#'
#' ## my.full
#' # $`1`
#' #     x1  x2  x3
#' # x1 1.0 0.3 0.4
#' # x2 0.3 1.0 0.5
#' # x3 0.4 0.5 1.0
#' #
#' # $`2`
#' #     x1 x2  x3
#' # x1 1.0 NA 0.4
#' # x2  NA NA  NA
#' # x3 0.4 NA 1.0
#'
#' ## Write two lower triangle correlation matrices into a file named "lowertriangle.dat".
#' ## x2 is missing in the second matrix.
#' ## The content of "lowertriangle.dat" is
#' # 1.0 
#' # 0.3 1.0 
#' # 0.4 0.5 1.0
#' # 1.0
#' # NA NA 
#' # 0.4 NA 1.0
#' ## cat("1.0\n0.3 1.0\n0.4 0.5 1.0\n1.0\nNA NA\n0.4 NA 1.0",
#' ##     file="lowertriangle.dat", sep="")
#'
#' ## Read the lower triangle correlation matrices from a file
#' ## my.lowertri <- readLowTriMat(file = "lowertriangle.dat", no.var = 3)
#'
#' ## Read the correlation matrices from a string
#' x <-
#' "1.0 
#' 0.3 1.0 
#' 0.4 0.5 1.0
#' 1.0
#' NA NA 
#' 0.4 NA 1.0"
#'
#' my.lowertri <- readLowTriMat(textConnection(x), no.var = 3)
#'
#' ## my.lowertri
#' # $`1`
#' #     x1  x2  x3
#' # x1 1.0 0.3 0.4
#' # x2 0.3 1.0 0.5
#' # x3 0.4 0.5 1.0
#' #
#' # $`2`
#' #     x1 x2  x3
#' # x1 1.0 NA 0.4
#' # x2  NA NA  NA
#' # x3 0.4 NA 1.0
#'
#' ## Write two vectors of correlation coefficients based on
#' ##  column major into a file named "stackvec.dat".
#' ## x2 is missing in the second matrix.
#' ## The content of "stackvec.dat" is
#' # 1.0 0.3 0.4 1.0 0.5 1.0
#' # 1.0 NA 0.4 NA NA 1.0
#' ## cat("1.0 0.3 0.4 1.0 0.5 1.0\n1.0 NA 0.4 NA NA 1.0\n",
#' ##     file="stackvec.dat", sep="")
#'
#' ## Read the stack vectors from a file
#' ## my.vec <- readStackVec("stackvec.dat")
#'
#' ## Read the stack vectors from a string
#' x <- "
#' 1.0 0.3 0.4 1.0 0.5 1.0
#' 1.0 NA 0.4 NA NA 1.0"
#'
#' my.vec <- readStackVec(textConnection(x))
#'
#' ## my.vec
#' # $`1`
#' #     x1  x2  x3
#' # x1 1.0 0.3 0.4
#' # x2 0.3 1.0 0.5
#' # x3 0.4 0.5 1.0
#' #
#' # $`2`
#' #    x1 x2  x3
#' # x1 1.0 NA 0.4
#' # x2  NA NA  NA
#' # x3 0.4 NA 1.0
#' }
#'
NULL





#' Studies on Students' School Engagement and Achievement Reported by Roorda et
#' al. (2011)
#'
#' The data set includes 45 studies on the influence of affective
#' teacher-student relationships on students' school engagement and achievement
#' reported by Roorda et al. (2011).
#'
#' The variables are: \describe{ \item{data}{A list of 45 studies of
#' correlation matrices. The variables are \emph{pos} (positive teacher-student
#' relations), \emph{neg} (negative teacher-student relations), \emph{enga}
#' (student engagement), and \emph{achiev} (student achievement).} \item{n}{A
#' vector of sample sizes} \item{SES}{A vector of average socio-economic status
#' (SES) of the samples} }
#'
#' @name Roorda11
#' @docType data
#' @references Jak, S., & Cheung, M. W.-L. (2018). Addressing heterogeneity in
#' meta-analytic structural equation modeling using subgroup analysis.
#' \emph{Behavior Research Methods}, \bold{50}, 1359-1373.
#' @source Roorda, D. L., Koomen, H. M. Y., Spilt, J. L., & Oort, F. J. (2011).
#' The influence of affective teacher-student relationships on students' school
#' engagement and achievement a meta-analytic approach. \emph{Review of
#' Educational Research}, \emph{81}(4), 493-529.
#' @keywords datasets
#' @examples
#'
#' \donttest{
#'
#' ## Random-effects model: First stage analysis
#' random1 <- tssem1(Cov = Roorda11$data, n = Roorda11$n, method = "REM",
#'                   RE.type = "Diag")
#' summary(random1)
#'
#' ## Model specification in lavaan syntax
#' model <- "## Regression paths
#'           enga ~ b31*pos + b32*neg
#'           achiev ~ b43*enga
#'           ## Variances of pos and neg are fixed at 1
#'           pos ~~ 1*pos
#'           neg ~~ 1*neg
#'           ## Correlation between pos and neg
#'           pos ~~ p21*neg
#'           ## Error variances
#'           enga ~~ p33*enga
#'           achiev ~~ p44*achiev"
#'
#' RAM <- lavaan2RAM(model, obs.variables=c("pos", "neg", "enga", "achiev"))
#' RAM
#'
#' ## ## Equivalent RAM specification using create.mxMatrix()
#' ## varnames <- c("pos", "neg", "enga", "achiev")
#' ## A <- create.mxMatrix(c(0,0,0,0,
#' ##                        0,0,0,0,
#' ##                        "0.1*b31","0.1*b32",0,0,
#' ##                        0,0,"0.1*b43",0),
#' ##                      type = "Full", nrow = 4, ncol = 4, byrow = TRUE,
#' ##                      name = "A", as.mxMatrix = FALSE)
#' ## dimnames(A) <- list(varnames, varnames)
#' ## A
#' ##
#' ## S <- create.mxMatrix(c(1,
#' ##                        ".5*p21",1,
#' ##                        0,0,"0.6*p33",
#' ##                        0,0,0,"0.6*p44"),
#' ##                      type="Symm", byrow = TRUE,
#' ##                      name="S", as.mxMatrix = FALSE)
#' ## dimnames(S) <- list(varnames, varnames)
#' ## S
#'
#' ## Random-effects model: Second stage analysis
#' random2 <- tssem2(random1, RAM=RAM, diag.constraints=TRUE,
#'                   intervals="LB")
#' summary(random2)
#'
#' ## Display the model with the parameter estimates    
#' plot(random2)
#' }
#'
NULL





#' Correlation Matrices from Scalco et al. (2017)
#'
#' The data set includes correlation matrices using the theory of planned
#' behavior to predict organic food consumption reported by Scalco17 et al.
#' (2017).
#'
#' A list of data with the following structure: \describe{ \item{data}{A list
#' of correlation matrices. The variables are \emph{ATT} (attitude), \emph{SN}
#' (subjective norm), \emph{PBC} (perceived behavior control), \emph{BI}
#' (behavioral intention), and \emph{BEH} (behavior)} \item{n}{A vector of
#' sample sizes} \item{Age}{A vector of the mean age of the samples}
#' \item{Female}{A vector of the percentage of the female samples} }
#'
#' @name Scalco17
#' @docType data
#' @source Scalco, A., Noventa, S., Sartori, R., & Ceschi, A. (2017).
#' Predicting organic food consumption: A meta-analytic structural equation
#' model based on the theory of planned behavior. \emph{Appetite}, \bold{112},
#' 235-248.
#' @keywords datasets
#' @examples
#'
#' data(Scalco17)
#'
NULL





#' Correlations from Stadler et al. (2015)
#'
#' The data set includes correlations between complex problem solving and
#' intelligence reported by Stadler et al. (2015).
#'
#' A list of data with the following structure: \describe{ \item{ID}{ID of the
#' effect sizes} \item{Authors}{Authors of the studies} \item{Year}{Year of the
#' studies} \item{N}{Sample size} \item{CPSMeasure}{Complex problem solving
#' (CPS) measure} \item{IntelligenceMeasure}{Intelligence measure}
#' \item{r}{Correlation between CPS and intelligence} \item{v}{Sampling
#' variance of r} }
#'
#' @name Stadler15
#' @docType data
#' @source Stadler, M., Becker, N., Godker, M., Leutner, D., & Greiff, S.
#' (2015). Complex problem solving and intelligence: A meta-analysis.
#' \emph{Intelligence}, \bold{53}, 92-101.
#' @keywords datasets
NULL





#' Summary Method for tssem1, wls, meta, and meta3LFIML Objects
#'
#' It summaries results for various class.
#'
#'
#' @name summary
#' @aliases summary.tssem1FEM summary.tssem1FEM.cluster summary.tssem1REM summary.wls summary.wls.cluster summary.meta summary.meta3LFIML summary.reml summary.CorPop summary.Cor3L summary.bootuniR2 summary.osmasem summary.osmasem2 summary.mxsem print.summary.tssem1FEM print.summary.tssem1FEM.cluster print.summary.wls print.summary.meta print.summary.meta3LFIML print.summary.reml print.summary.CorPop print.summary.Cor3L print.summary.bootuniR2 print.summary.mxsem
#' @param object An object returned from either class \code{tssem1FEM}, class
#' \code{tssem1FEM.cluster}, class \code{tssem1REM}, class \code{wls}, class
#' \code{wls.cluster}, class \code{meta}, class \code{meta3LFIML}, class
#' \code{reml}, class \code{mxsem} or class \code{CorPop}.
#' @param x An object returned from either class \code{summary.tssem1FEM},
#' class \code{tssem1FEM.cluster}, class \code{summary.wls}, class
#' \code{summary.meta}, class \code{summary.meta3LFIML}, class
#' \code{summary.reml} or class \code{summary.CorPop}.
#' @param homoStat Logical. Whether to conduct a homogeneity test on the effect
#' sizes.
#' @param allX Logical. Whether to report the predictors and the auxiliary
#' variables.
#' @param robust Logicial. Whether to use robust standard error from
#' \code{\link[OpenMx]{imxRobustSE}}.
#' @param df.adjustment Numeric. Adjust the degrees of freedom manually. It may
#' be necessary if the df calculated is incorrect when
#' \code{diag.constraints=TRUE}.
#' @param probs Quantiles for the parameter estimates.
#' @param cutoff.chisq.pvalue Cutoff of the p-value for the chi-square
#' statistic.
#' @param cutoff.CFI The cutoff of the CFI.
#' @param cutoff.SRMR The cutoff of the SRMR.
#' @param cutoff.RMSEA The cutoff of the RMSEA.
#' @param fitIndices Whether to calculate the chi-square statistic and various
#' goodness-of-fit indices in osmasem. Note. It may take a while since
#' statistics of the saturated and independence models are required.
#' @param numObs The number of observations in calculating the fit statistics
#' in osmasem. If it is missing, the total number of observations is used.
#' @param \dots Further arguments to be passed to
#' \code{\link[stats]{printCoefmat}}
#' @note If the OpenMx status1 is either 0 or 1, the estimation is considered
#' fine. If the OpenMx status1 is other values, it indicates estimation
#' problems. Users should refer to `OpenMx` website for more details.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
#' \code{\link[metaSEM]{meta}}, \code{\link[metaSEM]{reml}},
#' \code{\link[metaSEM]{rCor}}, \code{\link[metaSEM]{bootuniR2}},
#' \code{\link[metaSEM]{osmasem}}
#' @keywords methods
NULL





#' Correlation coefficients reported by Tenenbaum and Leaper (2002)
#'
#' Forty-eight studies reported by Tenenbaum and Leaper (2002, Table 1).
#'
#' The variables are: \describe{ \item{Authors}{Authors of the study}
#' \item{Year}{Year of publication} \item{N}{Sample size} \item{r}{Correlation
#' between parents' gender schemas and their offspring's gender-related
#' cognitions.} \item{v}{Sampling variance of r}
#' \item{Publication_source}{Publication source: 1="top-tier journal",
#' 2="second-tier journal or book chapter", 3="dissertation", 4="other
#' unpublished study"} \item{Author_gender}{Gender of the first author:
#' "W"="woman", "M"="man"} \item{Parent_type}{Parent type: "M"="mother",
#' "F"="father", "MF"="mother and father"} \item{Parent_predictor}{Parent
#' predictor: "S"="self gender schema", "A"="gender attitudes about others"}
#' \item{Offspring_age}{Offspring age (months)} \item{Offspring_type}{Offspring
#' type: "D"="daughter", "S"="son", "DS"="daughter and son"}
#' \item{Offspring_outcome}{Offspring outcome: "S"="gender schema for self",
#' "A"="gender attitudes toward others", "I"="gender-related interests and
#' preferences", "W"="work-related attitudes"} }
#'
#' @name Tenenbaum02
#' @docType data
#' @source Tenenbaum, H. R., & Leaper, C. (2002). Are parents' gender schemas
#' related to their children's gender-related cognitions? A meta-analysis.
#' \emph{Developmental Psychology}, \emph{38}(4), 615-630.
#' https://doi.org/10.1037/0012-1649.38.4.615
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(Tenenbaum02)
#' }
#'
NULL





#' Second Stage analysis of the univariate R (uniR) approach
#'
#' It conducts the second stage analysis of the uniR analysis by fitting
#' structural equation models on the average correlation matrix.
#'
#' This function implements the univariate r approach proposed by Viswesvaran
#' and Ones (1995) to conduct meta-analytic structural equation modeling
#' (MASEM). It treats the average correlation matrix as if it was a covariance
#' matrix in fitting structural equation models. The harmonic mean of the
#' sample sizes in combining correlation coefficients is used as the sample
#' size in fitting structural equation models. It is included in this package
#' for research interests. The two-stage structural equation modeling (TSSEM)
#' approach is preferred (e.g., Cheung, 2015; Cheung & Chan, 2005).
#'
#' @name uniR2
#' @aliases uniR2mx uniR2lavaan
#' @param x An object of class \code{uniR1} from \code{\link[metaSEM]{uniR1}}.
#' @param RAM A RAM object including a list of matrices of the model returned
#' from \code{\link[metaSEM]{lavaan2RAM}}.
#' @param Amatrix If \code{RAM} is not specified, an \code{Amatrix} is
#' required. An asymmetric matrix in the RAM specification with
#' \code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be
#' converted into \code{\link[OpenMx]{MxMatrix-class}} by the
#' \code{as.mxMatrix} function.
#' @param Smatrix If \code{RAM} is not specified, an \code{Smatrix} is
#' required. A symmetric matrix in the RAM specification with
#' \code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be
#' converted into \code{\link[OpenMx]{MxMatrix-class}} by the
#' \code{as.mxMatrix} function.
#' @param Fmatrix If \code{RAM} is not specified, an \code{Fmatrix} is
#' required. A filter matrix in the RAM specification with
#' \code{\link[OpenMx]{MxMatrix-class}}. If it is \code{NULL} (the default), an
#' identity matrix with the same dimensions of \code{Cov} will be created. If
#' it is a matrix, it will be converted into
#' \code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix} function. It
#' is not required when there is no latent variable.
#' @param model.name A string for the model name in
#' \code{\link[OpenMx]{mxModel}}. If it is missing, the default is "UniR2".
#' @param suppressWarnings Logical. If \code{TRUE}, warnings are suppressed. It
#' is passed to \code{\link[OpenMx]{mxRun}}.
#' @param silent Logical. An argument to be passed to
#' \code{\link[OpenMx]{mxRun}}
#' @param run Logical. If \code{FALSE}, only return the mx model without
#' running the analysis.
#' @param model A model specified using lavaan syntax see
#' \code{\link[lavaan]{model.syntax}}
#' @param \dots Further arguments to be passed to either
#' \code{\link[OpenMx]{mxRun}} or \code{\link[lavaan]{sem}}. For
#' \code{\link[lavaan]{sem}}, \code{fixed.x=FALSE} is passed automatically.
#' @return A fitted object returned from \code{\link[OpenMx]{mxRun}} or
#' \code{\link[lavaan]{sem}}.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{uniR1}}, \code{\link[metaSEM]{lavaan2RAM}},
#' \code{\link[metaSEM]{Becker09}}
#' @references Cheung, M. W.-L. (2015). \emph{Meta-analysis: A structural
#' equation modeling approach}. Chichester, West Sussex: John Wiley & Sons,
#' Inc.
#'
#' Cheung, M. W.-L., & Chan, W. (2005). Meta-analytic structural equation
#' modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10},
#' 40-64.
#'
#' Viswesvaran, C., & Ones, D. S. (1995). Theory testing: Combining
#' psychometric meta-analysis and structural equations modeling.
#' \emph{Personnel Psychology}, \bold{48}, 865-885.
#' @keywords uniR
NULL





#' Dataset on the effectiveness of multidimensional family therapy in treating
#' adolescents with multiple behavior problems
#'
#' This dataset includes 61 effect sizes from 19 manuscripts nested from 8
#' studies reported by van der Pol et al. (2017). It studies the effectiveness
#' of multidimensional family therapy in treating adolescents with multiple
#' behavior problems.
#'
#' A list of data with the following structure: \describe{ \item{Number}{Number
#' of the effect size.} \item{Study}{Authors of the studies.} \item{N}{Total
#' sample size.} \item{N_target}{Sample size in the target group.}
#' \item{N_control}{Sample size in the control group.}
#' \item{Comparison_condition}{Either cognitive behavioral therapy
#' (\code{CBT}), combined treatment (\code{CT}) or group therapy
#' (\code{Group}).} \item{Study_ID}{Level-3 cluster.} \item{Age_mean}{Mean age
#' of the participants.} \item{Fllow_up}{Follow-up duration (in months).}
#' \item{Per_Males}{Percentage of males.} \item{Per_Minorities}{Percentage of
#' minorities.} \item{Per_Conduct_disorder}{Percentage of participants with
#' conduct disorder} \item{Per_Severe_cannabis_users}{Percentage of
#' participants of severe cannabis use.} \item{Outcome_measure}{Either
#' substance abuse, delinquency, externalizing and internalizing
#' psychopathology, and family functioning} \item{d}{Effect size in Cohen's d.}
#' \item{v}{Sampling variance of d} }
#'
#' @name vanderPol17
#' @docType data
#' @source van der Pol, T. M., Hoeve, M., Noom, M. J., Stams, G. J. J. M.,
#' Doreleijers, T. A. H., van Domburgh, L., & Vermeiren, R. R. J. M. (2017).
#' Research Review: The effectiveness of multidimensional family therapy in
#' treating adolescents with multiple behavior problems - a meta-analysis.
#' \emph{Journal of Child Psychology and Psychiatry}, \bold{58}(5), 532-545.
#' https://doi.org/10.1111/jcpp.12685
#' @keywords datasets
#' @examples
#'
#' data(vanderPol17)
#'
NULL





#' Extract Covariance Matrix Parameter Estimates from Objects of Various
#' Classes
#'
#' It extracts the variance-covariance matrix of the parameter estimates from
#' objects of various classes.
#'
#'
#' @name vcov
#' @aliases vcov.tssem1FEM vcov.tssem1FEM.cluster vcov.tssem1REM vcov.wls vcov.wls.cluster vcov.meta vcov.meta3LFIML vcov.reml vcov.osmasem vcov.osmasem2 vcov.mxsem
#' @param object An object returned from objects of various classes
#' @param select Select \code{all} for both fixed- and random-effects
#' parameters, \code{fixed} for the fixed-effects parameters or \code{random}
#' for the random-effects parameters. For \code{meta3LFIML} objects,
#' \code{allX} is used to extract all parameters including the predictors and
#' auxiliary variables.
#' @param robust Logicial. Whether to use robust standard error from
#' \code{\link[OpenMx]{imxRobustSE}}.
#' @param \dots Further arguments; currently not in use except for
#' \code{tssemRobust1}, which to be passed to \code{\link[metafor]{robust}}.
#' @return A variance-covariance matrix of the parameter estimates.
#' @note \code{vcov} returns \code{NA} when the \code{diag.constraints=TRUE}
#' argument is used in \code{wls} objects.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
#' \code{\link[metaSEM]{meta}}, \code{\link[metaSEM]{reml}}
#' @keywords methods
#' @examples
#'
#' ## Random-effects meta-analysis
#' model1 <- meta(y=yi, v=vi, data=Hox02)
#' vcov(model1)
#'
#' ## Fixed-effects only
#' vcov(model1, select="fixed")
#'
#' ## Random-effects only
#' vcov(model1, select="random")
#'
NULL





#' Forty-four Studies from Cheung (2013)
#'
#' Between 1990 and 1993, 57,561 adults aged 18 and above from 42 nations were
#' interviewed by local academic institutes in Eastern European nations and by
#' professional survey organizations in other nations. The standardized mean
#' difference (SMD) between males and females on life satisfaction and life
#' control in each country were calculated as the effect sizes. Positive values
#' indicate that males have higher scores than females do.
#'
#' The variables are: \describe{ \item{country}{Country} \item{lifesat}{SMD on
#' life satisfaction} \item{lifecon}{SMD on life control}
#' \item{lifesat_var}{Sampling variance of lifesat} \item{inter_cov}{Sampling
#' covariance between lifesat and lifecon} \item{lifecon_var}{Sampling variance
#' of lifecon} \item{gnp}{Gross National Product} }
#'
#' @name wvs94a
#' @docType data
#' @references Au, K., & Cheung, M. W.-L. (2004). Intra-cultural variation and
#' job autonomy in 42 countries. \emph{Organization Studies}, \bold{25},
#' 1339-1362.
#'
#' Cheung, M. W.-L. (2013). Multivariate meta-analysis as structural equation
#' models. \emph{Structural Equation Modeling}, \bold{20}, 429-454.
#' @source World Values Study Group. (1994). World Values Survey, 1981-1984 and
#' 1990-1993 [Computer file]. \emph{Ann Arbor, MI: Inter-university Consortium
#' for Political and Social Research.}
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(wvs94a)
#'
#' ## Random-effects model
#' random.ma1 <- meta(y=cbind(lifesat, lifecon),
#'                    v=cbind(lifesat_var, inter_cov, lifecon_var), data=wvs94a,
#'                    model.name="Random effects model")
#' summary(random.ma1)
#'
#' ## Random-effects model with both population effect sizes fixed at 0
#' random.ma2 <- meta(y=cbind(lifesat, lifecon),
#'                    v=cbind(lifesat_var, inter_cov, lifecon_var), data=wvs94a,
#'                    intercept.constraints=matrix(0, nrow=1, ncol=2),
#'                    model.name="Effect sizes are fixed at 0")
#' summary(random.ma2)
#'
#' ## Compare the nested models
#' anova(random.ma1, random.ma2)
#'
#' ## Fixed-effects model by fixing the variance component at 0 
#' fixed.ma <- meta(y=cbind(lifesat, lifecon),
#'                  v=cbind(lifesat_var, inter_cov, lifecon_var), data=wvs94a,
#'                  RE.constraints=matrix(0, ncol=2, nrow=2),
#'                  model.name="Fixed effects model")
#' summary(fixed.ma)
#'
#' ## Mixed-effects model
#' ## gnp is divided by 10000 and centered by using 
#' ## scale(gnp/10000, scale=FALSE)
#' mixed.ma1 <- meta(y=cbind(lifesat, lifecon),
#'                   v=cbind(lifesat_var, inter_cov, lifecon_var),
#'                   x=scale(gnp/10000, scale=FALSE), data=wvs94a,
#'                   model.name="GNP as a predictor")
#' summary(mixed.ma1)
#'
#' ## Mixed-effects model with equal regression coefficients
#' mixed.ma2 <- meta(y=cbind(lifesat, lifecon),
#'                   v=cbind(lifesat_var, inter_cov, lifecon_var),
#'                   x=scale(gnp/10000, scale=FALSE), data=wvs94a,
#'                   coef.constraints=matrix(c("0.0*Eq_slope",
#'                                             "0.0*Eq_slope"), nrow=2),
#'                   model.name="GNP as a predictor with equal slope")
#' summary(mixed.ma2)
#'
#' ## Compare the nested models
#' anova(mixed.ma1, mixed.ma2)
#'
#' ## Plot the multivariate effect sizes
#' plot(random.ma1, main="Estimated effect sizes and their 95% confidence ellipses",
#'      axis.label=c("Gender difference on life satisfaction",
#'                   "Gender difference on life control"))
#' }
#'
NULL





#' Forty-four Covariance Matrices on Life Satisfaction, Job Satisfaction, and
#' Job Autonomy
#'
#' Between 1990 and 1993, 57,561 adults aged 18 and above from 42 nations were
#' interviewed by local academic institutes in Eastern European nations and by
#' professional survey organizations in other nations. The covariance matrices
#' among Life Satisfaction, Job Satisfaction, and Job Autonomy were calculated.
#'
#' The variables are: \describe{ \item{data}{Covariance matrix among Life
#' Satisfaction (LS), Job Satisfaction (JS), and Job Autonomy (JA)}
#' \item{n}{Sample size in the country} }
#'
#' @name wvs94b
#' @docType data
#' @references Au, K., & Cheung, M. W.-L. (2004). Intra-cultural variation and
#' job autonomy in 42 countries. \emph{Organization Studies}, \bold{25},
#' 1339-1362.
#'
#' Cheung, M.W.-L., & Cheung, S.-F. (2016). Random-effects models for
#' meta-analytic structural equation modeling: Review, issues, and
#' illustrations. \emph{Research Synthesis Methods}, \bold{7}, 140-155.
#' @source World Values Study Group. (1994). World Values Survey, 1981-1984 and
#' 1990-1993 [Computer file]. \emph{Ann Arbor, MI: Inter-university Consortium
#' for Political and Social Research.}
#' @keywords datasets
#' @examples
#'
#' \donttest{
#' data(wvs94b)
#'
#' ## Get the indirect and the direct effects and
#' ## their sampling covariance matrices for each study
#' indirect1 <- indirectEffect(wvs94b$data, wvs94b$n)
#' indirect1
#'
#' ## Multivariate meta-analysis on the indirect and direct effects
#' indirect2 <- meta(indirect1[, c("ind_eff", "dir_eff")],
#'                   indirect1[, c("ind_var", "ind_dir_cov", "dir_var")])
#'
#' summary(indirect2)  
#' }
#'
NULL



