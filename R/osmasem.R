## Reorder the observed variables of RAM (the rows of F) to 'obslabels',
## the order of the variables in the data. The model-implied correlations
## are matched to the data by position, so a RAM in another order would pair
## them with the wrong coefficients. A RAM whose names share none with the
## data (e.g., A1, D2, ... vs. x1, x2, ...) is matched by position, as before
## 1.6.1, with a warning.
.match.RAM.obs <- function(RAM, obslabels, subset.hint=FALSE) {
  obs <- rownames(RAM$F)
  if (is.null(obs)) {
    if (nrow(RAM$F) != length(obslabels))
      stop("RAM has ", nrow(RAM$F), " observed variables, but the data have ",
           length(obslabels), ".\n")
    return(RAM)
  }
  if (length(obs) == length(obslabels) && !any(obs %in% obslabels)) {
    warning("None of the observed variables in RAM (", paste(obs, collapse=", "),
            ") are in the data (", paste(obslabels, collapse=", "),
            "), so they are matched by position (", obs[1], " with ",
            obslabels[1], ", ...). Rename the variables in RAM or in the data ",
            "to match them by name.", call.=FALSE)
    return(RAM)
  }
  if (length(obs) != length(obslabels) || !setequal(obs, obslabels))
    stop("The observed variables in RAM (", paste(obs, collapse=", "),
         ") do not match those in the data (", paste(obslabels, collapse=", "),
         "). ", if (subset.hint) "Use 'subset.variables' or a" else "Use a",
         " RAM with the same variables.\n")
  RAM$F <- RAM$F[obslabels, , drop=FALSE]
  RAM
}

## Check 'subset.variables' of osmasem() and osmasem2() against the
## variables in the data.
.check.subset.variables <- function(subset.variables, obslabels) {
  absent <- setdiff(subset.variables, obslabels)
  if (length(absent))
    stop(paste(absent, collapse=", "), " are not in the data.\n")
  duplicated.names <- unique(subset.variables[duplicated(subset.variables)])
  if (length(duplicated.names))
    stop("'subset.variables' contains ", paste(duplicated.names, collapse=", "),
         " more than once.\n")
  invisible(NULL)
}

## Warn if an estimated heterogeneity (random-effects) covariance matrix is
## not positive semidefinite. Bounds on its variances do not prevent this
## when the covariances are free: the likelihood requires only Tau+V to be
## positive definite.
.warn.Tau.not.psd <- function(Tau, name="Tau") {
  Tau <- as.matrix(Tau)
  if (nrow(Tau) < 2 || anyNA(Tau)) return(invisible(NULL))
  ev <- eigen((Tau + t(Tau))/2, symmetric=TRUE, only.values=TRUE)$values
  if (min(ev) < -1e-6 * max(1, abs(ev)))
    warning("The estimated heterogeneity matrix '", name, "' is not positive ",
            "semidefinite (smallest eigenvalue = ", signif(min(ev), 3), "), so ",
            "it is not a valid covariance matrix. Consider a simpler ",
            "random-effects structure, e.g., a diagonal one.", call.=FALSE)
  invisible(NULL)
}

## A function to extract levels of path directions started from IVs to DVs
.pathLevels <- function(A) {
  edges <- A$free | (!is.na(A$values) & A$values != 0)
  p <- ncol(edges)
  remaining <- rep(TRUE, p)
  levels <- list()
  while (any(remaining)) {
    eligible <- remaining & rowSums(edges[, remaining, drop=FALSE]) == 0
    if (!any(eligible))
      stop("The path matrix contains a directed cycle; implicit error variances require an acyclic model.")
    levels[[length(levels)+1L]] <- eligible
    remaining[eligible] <- FALSE
  }
  list(IV=levels[1L], DV=levels[-1L])
}



## A0: Intercept
## A1: 1st moderator
## Amatrix = A0 + A1*Ax[[1]] + A1*Ax[[2]]
## Smatrix = S0 + S1*Sx[[1]] + S2*Sx[[2]]


#' Create a model implied correlation matrix with implicit diagonal constraints
#' 
#' It creates implicit diagonal constraints on the model implied correlation
#' matrix by treating the error variances as functions of other parameters.
#' 
#' 
#' @param A0 A Amatrix, which will be converted into \code{MxMatrix-class} via
#' \code{as.mxMatrix}.
#' @param S0 A Smatrix, which will be converted into \code{MxMatrix-class} via
#' \code{as.mxMatrix}.
#' @param F0 A Fmatrix, which will be converted into \code{MxMatrix-class} via
#' \code{as.mxMatrix}. If omitted, an identity filter is used with the
#' variable names of \code{A0}, when available.
#' @param Ax A Amatrix of a list of Amatrix with definition variables as the
#' moderators of the Amatrix.
#' @param Sx A Smatrix of a list of Smatrix with definition variables as the
#' moderators of the Smatrix.
#' @param A.lbound A matrix of lower bound of the Amatrix. If a scalar is
#' given, the lbound matrix will be filled with this scalar.
#' @param A.ubound A matrix of upper bound of the Amatrix. If a scalar is
#' given, the ubound matrix will be filled with this scalar.
#' @return A list of \code{MxMatrix-class}. The model implied correlation
#' matrix is computed in \code{impliedR} and \code{vechsR}.
#' @note Since \code{A0} are the intercepts and \code{Ax} are the regression
#' coefficients, the parameters in \code{Ax} must be a subset of those in
#' \code{A0}. Path models must be acyclic when residual variances are derived
#' implicitly. Fixed dependent-variable residual variances in \code{S0}
#' are ignored with a warning.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{osmasem}}, \code{\link[metaSEM]{create.Tau2}},
#' \code{\link[metaSEM]{create.V}}
#' @keywords osmasem
#' @examples
#' 
#' \donttest{
#' ## Proposed model
#' model1 <- 'W2 ~ w2w*W1 + s2w*S1
#'            S2 ~ w2s*W1 + s2s*S1
#'            W1 ~~ w1WITHs1*S1
#'            W2 ~~ w2WITHs2*S2
#'            W1 ~~ 1*W1
#'            S1 ~~ 1*S1
#'            W2 ~~ Errw2*W2
#'            S2 ~~ Errs2*S2'
#' 
#' ## Convert into RAM    
#' RAM1 <- lavaan2RAM(model1, obs.variables=c("W1", "S1", "W2", "S2"))
#' 
#' ## No moderator    
#' M0 <- create.vechsR(A0=RAM1$A, S0=RAM1$S, F0=NULL, Ax=NULL, Sx=NULL)
#' 
#' ## Lag (definition variable) as a moderator on the paths in the Amatrix    
#' Ax <- matrix(c(0,0,0,0,
#'                0,0,0,0,
#'                "0*data.Lag","0*data.Lag",0,0,
#'                "0*data.Lag","0*data.Lag",0,0),
#'              nrow=4, ncol=4, byrow=TRUE)
#'                 
#' M1 <- create.vechsR(A0=RAM1$A, S0=RAM1$S, F0=NULL, Ax=Ax, Sx=NULL)    
#' 
#' ## Lag (definition variable) as a moderator on the correlation in the Smatrix
#' Sx <- matrix(c(0,"0*data.Lag",0,0,
#'                "0*data.Lag",0,0,0,
#'                0,0,0,"0*data.Lag",
#'                0,0,"0*data.Lag",0),
#'              nrow=4, ncol=4, byrow=TRUE)
#' 
#' M2 <- create.vechsR(A0=RAM1$A, S0=RAM1$S, F0=NULL, Ax=NULL, Sx=Sx)
#' }
#' 
create.vechsR <- function(A0, S0, F0=NULL, Ax=NULL, Sx=NULL,
                          A.lbound=NULL, A.ubound=NULL) {

  for (mat in c(list(A0, S0, F0), if (is.list(Ax)) Ax else list(Ax),
                if (is.list(Sx)) Sx else list(Sx))) .check.model.matrix(mat)

  if (is.matrix(A0)) {
    A0 <- as.mxMatrix(A0, name="A0")
  } else {
    A0@name <- "A0"
  }    

  if (is.matrix(S0)) {
    S0 <- as.mxMatrix(S0, name="S0")
  } else {
    S0@name <- "S0"
  }

  ## Basic checking
  checkRAM(Amatrix=A0, Smatrix=S0, cor.analysis=TRUE)
  ## ## Additional check on whether the diagonals of S0 are 0
  ## if (any(Diag(S0$free)==TRUE)) {
  ##     stop("Diagonal of 'S0' must be 0.\n")
  ## }       
    
  ## Create an Identity matrix
  Iden <- create.Fmatrix(rep(1, ncol(A0$values)), name="Iden")

  ## If F0 is not specified, assume all variables are selected.
  if (is.null(F0)) {
    Fmatrix <- Iden
    Fmatrix$name <- "Fmatrix"
    dimnames(Fmatrix@values) <- dimnames(Fmatrix@labels) <-
      dimnames(Fmatrix@free) <- dimnames(A0$values)
  } else {
    Fmatrix <- as.mxMatrix(F0, name="Fmatrix")
  }

  ## Consider error variances as functions of parameters
  ## Check types of variables
  ## dv (not ivs): variables pointed by some variables
  dv <- rowSums(A0$free | (!is.na(A0$values) & A0$values != 0)) > 0

  if (any(dv & !diag(S0$free)))
    warning("For correlation analysis, create.vechsR() computes dependent-variable residual variances; fixed values in S0 for those variables are ignored.", call.=FALSE)

  ## ## iv: variables pointed towards other variables
  ## iv <- apply(Amatrix$free, 2, any)
  ## ## med: mediators
  ## med <- iv & dv

  ## Error variances of dv are modified as 0.
  ## S0: Modified Smatrix without error variances on the dependent variables
  ## Fix the diagonal elements associated with dv to 0
  diag(S0$labels)[dv] <- NA
  diag(S0$values)[dv] <- 0
  diag(S0$free)[dv] <- FALSE

  ## ## Use it rather than Smatrix in mxCI()
  ## S1_labels <- S1$labels
  ## diag(S1_labels) <- NA
  ## S1_labels <- c(na.omit(unique(c(S1_labels))))

  ## Levels of directions, which variables are dv
  Level <- .pathLevels(A0)$DV

  ## When there is no predictor, only A0.
  ## Amatrix=A0
  A.mod <- Ax.mats <- S.mod <- Sx.mats <- list()
  if (is.null(Ax)) {
    Amatrix <- mxAlgebra(A0, name="Amatrix")
  } else {
    ## Convert Apred into a list if it is not
    if (!is.list(Ax)) Ax <- list(Ax)

    ## Basic checking on Ax
    lapply(Ax, function(x) checkRAM(Amatrix=x))

    for (i in seq_along(Ax)) {
      ## Parameters for the moderator based on the Ax, not A0
      ## When the labels include "data", they are labels of the moderator.
      index_mod <- apply(Ax[[i]], 1:2, function(x) grepl("data", x))

      ## A copy of A0 whose free parameters are based on Ax, not A0
      Ai <- A0
      Ai$free <- index_mod
      Ai$labels[!index_mod] <- NA
      Ai$name <- paste0("A", i)
      ## Labels = Labels in A0 + "_1", "_2", etc.
      Ai$labels <- apply(Ai$labels, 1:2, function(x) if (is.na(x)) NA else paste0(x, "_", i))
      A.mod[[i]] <- Ai

      ## Values of the moderators with definition variables
      Ax.mats[[i]] <- as.mxMatrix(Ax[[i]], name=paste0("Ax", i))
    }

    ## Amatrix= A0+A1*Ax1+A2*Ax2+...
    Amatrix <- mxAlgebraFromString(paste0("A0 + ", paste0("A", seq_along(Ax), "*Ax",
                                                         seq_along(Ax), collapse=" + ")),
                                   name="Amatrix")
  }

  ## When there is no predictor, only S0.
  ## SnoErr=S0
  if (is.null(Sx)) {
    SnoErr <- mxAlgebra(S0, name="SnoErr")
  } else {
    ## Convert Sx into a list if it is not
    if (!is.list(Sx)) Sx <- list(Sx)

    ## Basic checking on Sx
    ## Don't check the diagonals
    lapply(Sx, function(x) checkRAM(Smatrix=x, cor.analysis=FALSE))

    for (i in seq_along(Sx)) {
      ## Parameters for the moderator
      Si <- S0
      Si$name <- paste0("S", i)
      ## Labels = Labels in S0 + "_1", "_2", etc.
      Si$labels <- apply(Si$labels, 1:2, function(x) if (is.na(x)) NA else paste0(x, "_", i))
      ## Fix the diagonals of S1, S2... to 0
      values <- Si$values
      diag(values) <- 0
      Si$values <- values
      S.mod[[i]] <- Si

      ## Values of the moderators with definition variables
      Sx.mats[[i]] <- as.mxMatrix(Sx[[i]], name=paste0("Sx", i))
    }

    ## Not the final Smatrix as it does not include the error variances
    ## S_no_err= S0+S1*Sx1+S2*Sx2+...
    SnoErr <- mxAlgebraFromString(paste0("S0 + ", paste0("S", seq_along(Sx), "*Sx",
                                                        seq_along(Sx), collapse=" + ")),
                                  name="SnoErr")
  }

  ## No error variance involved for ALL variables
  invS0 <- mxAlgebra( solve(Iden-Amatrix)%&%SnoErr, name="invS0")

  ## Setup for error variances for mediators and dvs
  ## Excluded level1 as they are IVs
  F.list <- E.list <- invS.list <- list()
  for (i in seq_along(Level)) {
    ## Select the diagonals of mediators
    F.list[[i]] <- as.mxMatrix(diag(Level[[i]]), name=paste0("F", i))

    ## Extract the error variances of mediators based on the previous
    ## model implied covariance matrix (i-1)
    E.list[[i]] <- mxAlgebraFromString(paste0("vec2diag(diag2vec(Iden-invS", i-1,
                                              "))*F", i), name=paste0("E", i))

    ## Implied covariance matrix with estimated error variances on mediators for (i)
    invS.list[[i]] <- mxAlgebraFromString(paste0("solve(Iden-Amatrix)%&%(SnoErr+",
                                                 paste0("E", seq_len(i), collapse="+"), ")"),
                                          name=paste0("invS", i))
  }

  ## Final Smatrix=SnoErr +E1+E2+E3...
  ## Without dependent variables (e.g., a model of correlations only), there
  ## are no error variances to add.
  Smatrix <- mxAlgebraFromString(paste0("SnoErr", if (length(Level))
                                          paste0("+E", seq_along(Level), collapse="")),
                                 name="Smatrix")

  ## Final impliedR
  impliedR <- mxAlgebraFromString("(Fmatrix%*%solve(Iden-Amatrix))%&%Smatrix",
                                  name="impliedR")
  ## Vectorized form of impliedR
  vechsR <- mxAlgebra(t(vechs(impliedR)), name="vechsR")

  ## Output the matrices, followed by those of the moderators if there are any
  out <- c(list(Amatrix, Smatrix, Fmatrix, impliedR, vechsR, A0, S0, SnoErr, invS0, Iden),
           invS.list, E.list, F.list, A.mod, Ax.mats, S.mod, Sx.mats)

  ## Add names of the elements in the list
  names(out) <- sapply(out, function(x) x$name)

  ## Add lbound and ubound to A0
  if (!is.null(A.lbound)) {
    if (is.matrix(A.lbound)) {
      if (all(dim(A.lbound)==dim(out$A0$values))) {
        out$A0$lbound <- A.lbound
      } else {
        stop("The dimensions of 'A.lbound' do not match those in 'A0'.\n")
      }
    } else {  # lbound is not matrix
      out$A0$lbound  <- matrix(A.lbound, nrow=nrow(out$A0$values),
                               ncol=ncol(out$A0$values),
                               dimnames=dimnames(out$A0$values))
    }
  }

  if (!is.null(A.ubound)) {
    if (is.matrix(A.ubound)) {
      if (all(dim(A.ubound)==dim(out$A0$values))) {
        out$A0$ubound <- A.ubound
      } else {
        stop("The dimensions of 'A.ubound' do not match those in 'A0'.\n")
      }
    } else {  # ubound is not matrix
      out$A0$ubound  <- matrix(A.ubound, nrow=nrow(out$A0$values),
                               ncol=ncol(out$A0$values),
                               dimnames=dimnames(out$A0$values))
    } 
  }      
                
  out
}



#' Create a variance component of the heterogeneity of the random effects
#' 
#' It creates variance component of the heterogeneity of the random effects by
#' decomposing the variance component into matrices of correlation and standard
#' deviations.
#' 
#' 
#' @param RAM The RAM model for testing. \code{no.var} is calculated from it.
#' @param no.var If \code{RAM} is missing, the user has to specify the
#' \code{no.var} argument. It represents the \code{no.var} by \code{no.var} of
#' the random effects).
#' @param Tau1.labels Parameter labels in \code{Tau1}. The default is
#' \code{Tau1_1}, \code{Tau1_2}, etc. Labels must not contain periods.
#' @param RE.type Either \code{"Diag"}, \code{"CSymm"}, \code{"HCSymm"},
#' \code{"Symm"}, \code{"Zero"} or \code{"User"}. If it is\code{"Diag"} (the
#' default if missing), a diagonal matrix is used for the random effects
#' meaning that the random effects are independent. \code{"CSymm"} and
#' \code{"HCSymm"} represent compound symmetry and heteroscedastic compound
#' symmetry, respectively. If it is \code{"Symm"}, a symmetric matrix is used
#' for the random effects on the covariances among the correlation (or
#' covariance) vectors. If it is \code{"Zero"}, a zero matrix is assumed on the
#' variance component of the random effects. If it is \code{"User"}, users have
#' to specify the \code{RE.User} argument.
#' @param level whether it is for single-level, between-, or within-level
#' analyses. The only difference are the names of the matrices.
#' @param RE.User It represents the \code{no.var} by \code{no.var} symmetric
#' matrix of \code{TRUE} or \code{FALSE} for the variance component. If the
#' elements are \code{FALSE}, they are fixed at 0.
#' @param Transform Either \code{"expLog"} or \code{"sqSD"}. If it is
#' \code{"expLog"}, the variances are estimated by applying a log and exp
#' transformation. If it is \code{"sqSD"}, the variances are estimated by
#' applying a square on the SD. The transformation may improve the estimation
#' when the heterogeneity is small or close to zero.
#' @param RE.startvalues Starting values for the variances.
#' @return A list of \code{MxMatrix-class}. The variance component is computed
#' in \code{Tau2}.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{osmasem}}, \code{\link[metaSEM]{create.V}},
#' \code{\link[metaSEM]{create.vechsR}}
#' @keywords osmasem osmasem3L
#' @examples
#' 
#' \donttest{
#' T0 <- create.Tau2(no.var=4, RE.type="Diag", Transform="expLog", RE.startvalues=0.05)
#' T0
#' 
#' T1 <- create.Tau2(no.var=4, Tau1.labels=c("a", "b", "c", "d"))
#' T1
#' }
#' 
create.Tau2 <- function(RAM, no.var, Tau1.labels=seq(no.var),
                        RE.type=c("Diag", "CSymm", "HCSymm", "Symm", "Zero", "User"),
                        level=c("single", "between", "within"),
                        RE.User=NULL, Transform=c("expLog", "sqSD"),
                        RE.startvalues=0.05) {

  if (!missing(RAM)) .check.ram.names(RAM)

  level <- match.arg(level)
  switch(level,
         single = level_suffix <- "",
         between= level_suffix <- "B",
         within = level_suffix <- "W")

  if (!missing(RAM)) {
    p <- sum(RAM$F)
    ## no. of correlation coefficients
    no.var <- p*(p-1)/2
  }
    
  ## It has not been tested yet. But it should work for no.var=1.
  if (no.var<1) stop("'no.var' must be at least 1.\n")
  
  if (length(Tau1.labels) != no.var) stop("Length of \"Tau1.labels\" is different from \"no.var\".\n")
  .check.no.period(Tau1.labels, "Random-effect parameter labels")
  Tau1.labels <- paste0(Tau1.labels, level_suffix)
    
  RE.type <- match.arg(RE.type)
  Transform <- match.arg(Transform)

  ## Repeat the starting values if there is only one value.
  if (length(RE.startvalues)==1) {
    RE.startvalues <- rep(RE.startvalues, no.var)
  }
  
  ## Convert the starting values into positive.
  if (any(RE.startvalues < 0)) {
    RE.startvalues <- abs(RE.startvalues) + 1e-10
    warning("Not all 'RE.startvalues' are positive.\n")
  }

  ## Set the starting values at 0 for RE.type="Zero"
  if (RE.type=="Zero") RE.startvalues <- rep(0, no.var)
  
  ## Convert the starting values into log or sqrt scale.
  switch(Transform,
         expLog = RE.startvalues <- log(RE.startvalues),
         sqSD = RE.startvalues <- sqrt(RE.startvalues))

  ## vecTau1: px1 vector of the log (or SD) Tau1
  ## Cor: pxp corr matrix
  ## Tau2: pxp final variance component of variances
  switch(RE.type,
         Symm = {              
           vecTau1 <- create.mxMatrix(paste0(RE.startvalues, "*Tau1_", Tau1.labels),
                                      ncol=1, nrow=no.var, name="vecTau1")
           Cor <- create.mxMatrix(vechs(outer(seq(no.var), seq(no.var),
                                              function(x,y) paste0("0*Cor_", x, "_", y, level_suffix))),
                                  type="Stand", ncol=no.var, nrow=no.var,
                                  lbound=-0.99, ubound=0.99, name=paste0("Cor", level_suffix))},
         Diag = {
           vecTau1 <- create.mxMatrix(paste0(RE.startvalues, "*Tau1_", Tau1.labels),
                                      ncol=1, nrow=no.var, name="vecTau1")
           ## Uncorrelated
           Cor <- as.mxMatrix(diag(no.var), name=paste0("Cor", level_suffix))},
         Zero = {
           vecTau1 <- create.mxMatrix(RE.startvalues, type="Full", ncol=1,
                                      nrow=no.var, name="vecTau1")
           Cor <- as.mxMatrix(diag(no.var), name=paste0("Cor", level_suffix))},
         CSymm = {
           ## Compound Symmetry: One variance label for all, one correlation label for all
           vecTau1 <- create.mxMatrix(rep(paste0(RE.startvalues[1], "*Tau1_", Tau1.labels[1]), no.var),
                                      ncol=1, nrow=no.var, name="vecTau1")
           # Create a correlation matrix where all off-diagonals share the same label
           cor_labels <- matrix(NA, nrow=no.var, ncol=no.var)
           cor_labels[lower.tri(cor_labels)] <- paste0("Cor_csymm_", level_suffix)
           cor_labels[upper.tri(cor_labels)] <- t(cor_labels)[upper.tri(cor_labels)]
           Cor <- mxMatrix(type="Stand", nrow=no.var, ncol=no.var, free=TRUE, values=0,
                           labels=cor_labels,
                           lbound=-0.99, ubound=0.99, name=paste0("Cor", level_suffix))},
         HCSymm = {
           ## Heteroscedastic Compound Symmetry: Unique variances, one shared correlation
           vecTau1 <- create.mxMatrix(paste0(RE.startvalues, "*Tau1_", Tau1.labels),
                                      ncol=1, nrow=no.var, name="vecTau1")
           # Create a correlation matrix where all off-diagonals share the same label
           cor_labels <- matrix(NA, nrow=no.var, ncol=no.var)
           cor_labels[lower.tri(cor_labels)] <- paste0("Cor_hcsymm_", level_suffix)
           cor_labels[upper.tri(cor_labels)] <- t(cor_labels)[upper.tri(cor_labels)]
           Cor <- mxMatrix(type="Stand", nrow=no.var, ncol=no.var, free=TRUE, values=0,
                           labels=cor_labels,
                           lbound=-0.99, ubound=0.99, name=paste0("Cor", level_suffix))},       
         User = { if (is.null(RE.User)) stop("'RE.User', the ", no.var, " by ", no.var,
                                             " symmetric matrix of TRUE or FALSE on the variance componment, must be specified.\n")
                                             ## Check the symmetry of RE.User
                                             if (!isSymmetric(RE.User)) stop("'RE.User' is not symmetric.\n")
                                             if (!all(dim(RE.User)==c(no.var, no.var))) stop("The dimensions of 'RE.User' are different 'no.var'.\n")
                                             ## Check if the covariances are free but the variances are fixed.                                
                                             if (no.var>1) {
                                               for (i in 1:(no.var-1)) {
                                                 for (j in (i+1):no.var) {
                                                   if (RE.User[i,j] & !all(RE.User[i,i], RE.User[j,j])) {
                                                     stop("RE.User[",j,",",i,"] is free but either RE.User[",i,",",
                                                          i,"] or RE.User[",j,",",j,"] is fixed.\n") }}}}
                                               
                                             vecTau1 <- paste0(RE.startvalues, "*Tau1_", Tau1.labels)
                                             ## Fixed a bug that RE.startvalues are not used when the variances are not free.
                                             ## vecTau1[diag(RE.User)==FALSE] <- 0                                
                                             switch(Transform,
                                                    expLog = vecTau1[diag(RE.User)==FALSE] <- log(0),
                                                    sqSD =   vecTau1[diag(RE.User)==FALSE] <- sqrt(0))
                                             vecTau1 <- create.mxMatrix(vecTau1, ncol=1, nrow=no.var, name="vecTau1")
                                             Cor <- outer(seq(no.var), seq(no.var),
                                                          function(x,y) paste0("0*Cor_", x, "_", y, level_suffix))
                                             Cor[RE.User==FALSE] <- 0
                                             Cor <- create.mxMatrix(vechs(Cor), type="Stand", ncol=no.var, nrow=no.var,
                                                                    lbound=-0.99, ubound=0.99, name=paste0("Cor", level_suffix))}
         )
  
  switch(Transform,
         expLog = Tau2 <- eval(parse(text=paste0("mxAlgebra( vec2diag(sqrt(exp(vecTau1))) %&% ",
                                                 "Cor", level_suffix,
                                                 ", name='Tau2", level_suffix, "')"))) ,
         sqSD = Tau2 <- eval(parse(text=paste0("mxAlgebra( vec2diag(vecTau1) %&% ",
                                               "Cor", level_suffix,
                                               ", name='Tau2", level_suffix, "')"))))

  ## Vectorized version of Tau2
  vechTau2 <- eval(parse(text=paste0("mxAlgebra(vech(Tau2", level_suffix,
                                     "), name='vechTau2", level_suffix, "')")))
    
  out <- list(Tau2, vecTau1, Cor, vechTau2)

  ## Add names of the elements in the list
  names(out) <- sapply(out, function(x) x$name)
  out
}    

## x: labels of variables by column major


#' Create a V-known matrix
#' 
#' It creates a V-known matrix of the sampling covariance matrix using
#' definition variables.
#' 
#' 
#' @param x A character vector of variable names of the sampling covariance
#' matrix.
#' @param type Either \code{"Symm"}, \code{"Diag"} or \code{"Full"}. Suppose
#' the number of variables is \eqn{p}, the numbers of variable names for
#' \code{"Symm"}, \code{"Diag"}, and \code{"Full"} are \eqn{p(p-1)/2
#' }{p(p-1)/2}, \eqn{p}{p}, and \eqn{p*p}{p*p}, respectively. The elements are
#' arranged in a column major.
#' @param as.mxMatrix Logical. Whether to convert the output into
#' \code{MxMatrix-class}.
#' @return A list of \code{MxMatrix-class}. The V-known sampling covariance
#' matrix is computed in \code{V}.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{osmasem}}, \code{\link[metaSEM]{create.Tau2}},
#' \code{\link[metaSEM]{create.vechsR}}
#' @keywords osmasem
#' @examples
#' 
#' \donttest{
#' my.df <- Cor2DataFrame(Nohe15A1)
#' 
#' ## Create known sampling variance covariance matrix
#' V0 <- create.V(my.df$vlabels)
#' V0
#' }
#' 
create.V <- function(x, type=c("Symm", "Diag", "Full"), as.mxMatrix=TRUE) {
  type <- match.arg(type)

  switch(type,
         Symm = out <- vec2symMat(x, diag=TRUE, byrow=FALSE),
         Diag = out <- Diag(x),
         Full = {
           no.var <- sqrt(length(x))
           if (no.var%%1 != 0) {
             stop("Length of 'x' (", length(x), ") is not a square.\n")
           } else {
             out <- matrix(x, nrow=no.var, ncol=no.var, byrow=FALSE)
           }
         })
           
  out <- apply(out, c(1, 2),
               function(x) { if (x=="0") x else paste0("0*data.", x) })
  
  if (as.mxMatrix) out <- as.mxMatrix(out, name="V")
  out
}
    


#' One-stage meta-analytic structural equation modeling
#' 
#' It fits MASEM with the one-stage MASEM (OSMASEM) approach.
#' Variable names and data column names must not contain periods. Rename them consistently, for example \code{x.1} to \code{x_1}, and check that the new names are unique. Study, group, and cluster values may contain periods.
#' 
#' \code{osmasem()} was implemented based on Jak and Cheung (2020) for
#' meta-analyzing correlation matrices. \code{osmasem2()} was a rewritten
#' version to handle correlation or covariance matrices with the means. There
#' are several major differences between them: 1. \code{osmasem()} allows the
#' use of \code{RAM} or (\code{Mmatrix} and \code{Tmatrix}), while
#' \code{osmasem2()} only accepts the \code{RAM} input. 2. \code{RE.type} is
#' used to specify the type of random effects on the correlations in
#' \code{osmasem()}. On the contrary, \code{osmasem2()} has three types of
#' random effects: correlations/covariances, means, and covariance between
#' correlations/covariance and means. 3. \code{osmasem()} reports the
#' transformed random effects in the parameter table. Users have to use
#' \code{VarCorr()} to obtain the heterogeneity matrix of the random effects.
#' In contrast, \code{osmasem2()} reports the heterogeneity matrix in the
#' parameter table. 4. \code{osmasem()} accepts either
#' \code{intervals.type="z"} or \code{intervals.type="LB"}, whereas
#' \code{osmasem2()} only uses \code{intervals.type="z"}. Thus, this argument
#' is removed in \code{osmasem2()}: a request for \code{LB} raises an error and \code{z} is deprecated. 5. Both functions accept linear and
#' nonlinear constraints and parameter functions (\code{:=}) in \code{RAM};
#' only \code{osmasem2()} can substitute equality constraints into the model
#' with \code{replace.constraints=TRUE}.
#'
#' The observed variables of \code{RAM} are matched to the data by name and
#' reordered to the order of the variables in the data. They must be the
#' variables in \code{subset.variables} if it is given, and all variables in
#' the data otherwise. If none of the names of the observed variables in
#' \code{RAM} are in the data (e.g., \code{A1}, \code{D2}, ... in \code{RAM}
#' and \code{x1}, \code{x2}, ... in the data) and their numbers are equal,
#' they are matched by position with a warning. Both functions warn if the estimated
#' heterogeneity matrix is not positive semidefinite, which bounds on its
#' variances cannot prevent when its covariances or correlations are free.
#' 
#' @aliases osmasem osmasem2
#' @param model.name A string for the model name in
#' \code{\link[OpenMx]{mxModel}}.
#' @param RAM A RAM object including a list of matrices of the model returned
#' from \code{\link[metaSEM]{lavaan2RAM}}. If it is given, \code{Mmatrix} and
#' \code{Tmatrix} arguments will be ignored.
#' @param Mmatrix A list of matrices of the model implied correlation matrix
#' created by the \code{create.vechsR}. It is only required when \code{RAM} is
#' null. A named filter matrix is matched to the data by name, retaining its
#' parameter settings; an unnamed filter is matched by position. Its number
#' of rows must match the number of observed variables in the data.
#' @param Tmatrix A list of matrices of the heterogeneity variance-covariance
#' matrix created by the \code{create.Tau2}. It is only required when
#' \code{RAM} is null.
#' @param Jmatrix The Jacobian matrix of the mean structure in mxMatrix. The
#' covariance structure is Jmatrix \%&\% Tau2 + Vi. If it is not given, an
#' identity matrix will be used.
#' @param Ax A Amatrix of a list of Amatrix with definition variables as the
#' moderators of the Amatrix. It is used to create the \code{Mmatrix}.
#' @param Sx A Smatrix of a list of Smatrix with definition variables as the
#' moderators of the Smatrix. It is used to create the \code{Mmatrix}.
#' @param A.lbound A matrix of lower bound of the Amatrix. If a scalar is
#' given, the lbound matrix will be filled with this scalar. Bounds from
#' inequalities in the lavaan syntax of \code{RAM} (see
#' \code{\link{lavaan2RAM}}) are applied to the parameters with those labels
#' in any matrix; on the Amatrix, \code{A.lbound} and \code{A.ubound} take
#' precedence over them when the corresponding entry is not \code{NA};
#' \code{NA} entries retain RAM bounds. Conflicting explicit bounds for
#' repeated labels are rejected. Bounds on parameters that \code{osmasem()} computes
#' rather than estimates (e.g., error variances of dependent variables) are
#' ignored with a warning.
#' @param A.ubound A matrix of upper bound of the Amatrix. If a scalar is
#' given, the ubound matrix will be filled with this scalar.
#' @param RE.type Type of the random effects. See
#' \code{\link[metaSEM]{create.Tau2}} for details.
#' @param RE.User A symmetric matrix of \code{TRUE} or \code{FALSE} for the
#' variance component when \code{RE.type="User"}. See
#' \code{\link[metaSEM]{create.Tau2}} for details.
#' @param data A list of data created by the \code{Cor2DataFrame}.
#' @param subset.variables A character vector of the observed variables
#' selected for the analysis, in any order. The analysis keeps the order of
#' the variables in the data. Each variable can be given only once. In
#' \code{osmasem2()}, it also selects the means, their sampling covariances,
#' and the sampling covariances between the correlations or covariances and
#' the means (\code{VyMlabels}, if present). Studies without any of the
#' selected correlations or covariances do not contribute to the fit.
#' @param subset.rows A logical vector of the same length as the number of rows
#' in the data to select which rows are used in the analysis.
#' @param intervals.type Either \code{z} (default if missing) or \code{LB}. If
#' it is \code{z}, it calculates the 95\% confidence intervals (CIs) based on
#' the estimated standard error. If it is \code{LB}, it calculates the 95\%
#' likelihood-based CIs on the parameter estimates.
#' @param mxModel.Args A list of arguments passed to
#' \code{\link[OpenMx]{mxModel}}.
#' @param mxRun.Args A list of arguments passed to \code{\link[OpenMx]{mxRun}}.
#' @param suppressWarnings Logical. If it is \code{TRUE}, warnings are
#' suppressed. This argument is passed to \code{\link[OpenMx]{mxRun}}. If
#' it is \code{FALSE}, warnings during the fit are shown after it, so that
#' they cannot stop it. In both cases, they are stored in the \code{warning}
#' element of the output.
#' @param silent Logical. An argument is passed to \code{\link[OpenMx]{mxRun}}
#' @param run Logical. If \code{FALSE}, only return the mx model without
#' running the analysis.
#' @param \dots Not used, except that \code{osmasem()} rejects
#' \code{cor.analysis=FALSE}; use \code{osmasem2()} for covariance data.
#' @param cor.analysis Whether to analyze correlation or covariance structure
#' analysis.
#' @param RE.type.Sigma Type of the random effects of the correlation or
#' covariance vectors.
#' @param RE.type.Mu Type of the random effects of the mean vectors.
#' @param RE.type.SigmaMu Type of the random effects between the
#' correlation/covariance vectors and the mean vectors.
#' @param mean.analysis Whether to include the analysis of the mean structure.
#' @param startvalues An optional list of starting values. It is useful when
#' there are new parameters in RAM.
#' @param replace.constraints It is relevant only when there are constraints in
#' RAM. If it is \code{FALSE}, these constraints will be imposed. If it is
#' \code{TRUE}, the parameters on the left-hand side will be replaced by the
#' algebras on the right-hand side. In \code{osmasem2()}, chained replacements
#' are resolved before substitution, including in remaining constraints and
#' defined parameters; parameter bounds are retained.
#' Parameters appearing only in equality expressions require substitution or
#' an explicit parameter matrix supplied through \code{mxModel.Args}.
#' For correlation analysis in \code{osmasem2()}, unit-variance restrictions
#' are imposed for all selected observed variables at every distinct observed
#' definition-variable pattern, even when a residual variance is fixed.
#' Algebraic identities and linearly dependent polynomial restrictions are
#' omitted where they can be established for an acyclic RAM, allowing for
#' floating-point cancellation relative to each monomial's coefficients. Small
#' independent terms are retained. Residual variance
#' specifications must be compatible with the moderator effects across studies.
#' @return An object of class \code{osmasem}
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{Cor2DataFrame}},
#' \code{\link[metaSEM]{create.vechsR}}, \code{\link[metaSEM]{create.Tau2}},
#' \code{\link[metaSEM]{create.V}}, \code{\link[metaSEM]{osmasem}},
#' \code{\link[metaSEM]{Nohe15}}, \code{\link[metaSEM]{issp05}}
#' @references Jak, S., & Cheung, M. W.-L. (2020). Meta-analytic structural
#' equation modeling with moderating effects on SEM parameters.
#' \emph{Psychological Methods}, \bold{25} (4), 430-455.
#' https://doi.org/10.1037/met0000245
#' 
#' Jak, S., Cheung, M. W.-L., Acar, S., & Kindred, R. (December 18, 2024).
#' Evaluating differences in latent means across studies: Extending
#' meta-analytic confirmatory factor analysis with the analysis of means. OSF.
#' https://doi.org/10.31234/osf.io/35gtz.
#' @keywords meta-analysis
osmasem <- function(model.name="osmasem", RAM=NULL, Mmatrix=NULL,
                    Tmatrix=NULL, Jmatrix=NULL, Ax=NULL, Sx=NULL,
                    A.lbound=NULL, A.ubound=NULL,
                    RE.type=c("Diag", "CSymm", "HCSymm", "Symm", "Zero", "User"),
                    RE.User=NULL, data,
                    subset.variables=NULL, subset.rows=NULL, 
                    intervals.type = c("z", "LB"),
                    mxModel.Args=NULL, mxRun.Args=NULL,
                    suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {
  .check.ram.names(RAM)
  .check.no.period(data$obslabels, "Observed variable names")
  .check.no.period(colnames(data$data), "Data column names")
  .check.no.period(subset.variables, "subset.variables")
  if (!is.null(Mmatrix))
    for (mat in Mmatrix) .check.model.matrix(mat)
  dots <- list(...)
  if ("cor.analysis" %in% names(dots) && !identical(dots$cor.analysis, TRUE))
    stop("osmasem() supports correlation analysis only; use osmasem2(cor.analysis=FALSE) for covariances.", call.=FALSE)
  expected <- .genCorNames(data$obslabels, cor.analysis=TRUE)$ylabels
  if (length(data$ylabels) != length(expected) || anyDuplicated(data$ylabels) ||
      !setequal(data$ylabels, expected))
    stop("osmasem() requires correlation data from Cor2DataFrame(cor.analysis=TRUE).", call.=FALSE)

  RE.type <- match.arg(RE.type)
    
  intervals.type <- match.arg(intervals.type)
  switch(intervals.type,
         z = intervals <- FALSE,
         LB = intervals <- TRUE)

  ## Operator for not in
  ## "%ni%" <- Negate("%in%")

  ## Filter out variables not used in the analysis
  if (!is.null(subset.variables)) {
    .check.subset.variables(subset.variables, data$obslabels)
    ## Keep the order of the variables in the data, which the names of the
    ## correlations and their sampling covariances follow.
    temp <- .subset.Cor2DataFrame(data, subset.variables, cor.analysis=TRUE)
    obslabels <- temp$obslabels
    ylabels <- temp$ylabels
    vlabels <- temp$vlabels
  } else {
    obslabels <- data$obslabels
    ylabels <- data$ylabels
    vlabels <- data$vlabels
  }

  ## If RAM is provided, create the matrices based on it.
  if (!is.null(RAM)) {
    RAM <- .match.RAM.obs(RAM, obslabels, subset.hint=is.null(subset.variables))
    ## Bounds from inequalities in the lavaan syntax (see lavaan2RAM()) are
    ## applied by label to the model below. On the A matrix, 'A.lbound' and
    ## 'A.ubound' take precedence.
    A.parameters <- as.mxMatrix(RAM$A, name="A")
    A.labels <- A.parameters$labels
    RAM.lbound <- attr(RAM, "lbound")
    RAM.ubound <- attr(RAM, "ubound")
    overridden <- function(bound) {
      if (is.null(bound)) return(NULL)
      if (!(is.numeric(bound) || is.logical(bound)) || is.complex(bound) || !length(bound))
        stop("A bounds must be numeric and nonempty.", call.=FALSE)
      if (is.logical(bound)) bound <- bound * 1
      if (!is.matrix(bound))
        bound <- matrix(bound, nrow=nrow(A.labels), ncol=ncol(A.labels))
      if (!identical(dim(bound), dim(A.labels)))
        stop("The dimensions of the A bounds do not match those in the RAM A matrix.", call.=FALSE)
      use <- !is.na(bound) & !is.na(A.labels) & A.parameters$free
      labels <- A.labels[use]; values <- bound[use]
      for (lab in unique(labels))
        if (length(unique(values[labels == lab])) > 1L)
          stop("Conflicting A bounds for shared parameter '", lab, "'.", call.=FALSE)
      keep <- !duplicated(labels)
      as.list(setNames(values[keep], labels[keep]))
    }
    A.lower.override <- overridden(A.lbound)
    A.upper.override <- overridden(A.ubound)
    RAM.lbound <- RAM.lbound[setdiff(names(RAM.lbound), names(A.lower.override))]
    RAM.ubound <- RAM.ubound[setdiff(names(RAM.ubound), names(A.upper.override))]
    Mmatrix <- create.vechsR(A0=RAM$A, S0=RAM$S, F0=RAM$F, Ax=Ax, Sx=Sx,
                             A.lbound=A.lbound, A.ubound=A.ubound)
    Tmatrix <- create.Tau2(RAM=RAM, RE.type=RE.type, RE.User=RE.User, Transform="expLog")
  } else if (inherits(Mmatrix$Fmatrix, "MxMatrix")) {
    filter <- Mmatrix$Fmatrix
    if (nrow(filter$values) != length(obslabels))
      stop("Mmatrix has ", nrow(filter$values), " observed variables, but the data have ",
           length(obslabels), ". Check 'subset.variables' and the model.", call.=FALSE)
    obs <- rownames(filter$values)
    if (is.null(obs) && !is.null(Mmatrix$A0) &&
        nrow(filter$values) == ncol(filter$values) &&
        isTRUE(all.equal(unname(filter$values), diag(nrow(filter$values)))) &&
        !any(filter$free)) {
      obs <- rownames(Mmatrix$A0$values)
      if (!is.null(obs)) dimnames(filter@values) <- list(obs, obs)
    }
    if (!is.null(obs)) {
      if (length(obs) == length(obslabels) && !anyDuplicated(obs) &&
          !anyDuplicated(obslabels) && setequal(obs, obslabels)) {
        filter <- filter[match(obslabels, obs), , drop=FALSE]
      } else {
        warning("Observed variable names in Mmatrix and data do not match; ",
                "matching by position. Check the variable order and possible typos; ",
                "estimates may be affected.", call.=FALSE)
      }
    }
    Mmatrix$Fmatrix <- filter
  }

  ## Create known sampling variance covariance matrix
  Vmatrix <- create.V(vlabels, type="Symm", as.mxMatrix=TRUE)

  ## If Jmatrix is not specified, create an identify matrix
  ## Use Vmatrix$Cor$values as it is an identity matrix
  if (is.null(Jmatrix)) {
    Jmatrix <- Vmatrix$values
    diag(Jmatrix) <- 1
    Jmatrix <- as.mxMatrix(Jmatrix, name="Jmatrix")
  } else {
    Jmatrix$name  <- "Jmatrix"
  }
    
  if (is.null(subset.rows)) {
    subset.rows <- rep(TRUE, nrow(data$data))
  }
  ## Select data for the analaysis
  mx.data <- data$data[subset.rows, ]

  ## Dirty trick to avoid the warning of "no visible binding for global variable"
  mx.model <- eval(parse(text="mxModel(model=model.name,
                           mxData(observed=mx.data, type='raw'),
                           mxExpectationNormal(covariance='expCov',
                                               means='vechsR',
                                               dimnames=ylabels),
                           mxFitFunctionML(),
                           Mmatrix, Tmatrix, Vmatrix, Jmatrix,
                           mxAlgebra(Jmatrix %&% Tau2 + V, name='expCov'),
                           mxCI(c('Amatrix', 'Smatrix', 'Tau2')))"))
    
  ## Add additiona arguments from RAM to mxModel
  if (!is.null(RAM$mxalgebras)) {
    for (i in seq_along(RAM$mxalgebras)) {
      mx.model <- mxModel(mx.model, RAM$mxalgebras[[i]])
    }
    ## check if they are mxalgebra, not mxconstraint
    algebra.names <- names(RAM$mxalgebras)
    isalgebra <- !vapply(RAM$mxalgebras, inherits, logical(1),
                         what="MxConstraint")
    if (any(isalgebra)) {
      mx.model <- mxModel(mx.model, mxCI(algebra.names[isalgebra]))
    }
  }  

  ## Add additional arguments to mxModel
  if (!is.null(mxModel.Args)) {
    for (i in seq_along(mxModel.Args)) {
      mx.model <- mxModel(mx.model, mxModel.Args[[i]])
    }
  }

  ## Assign same starting values, e.g., A and S matrices are different from those in provided in Jmatrix
  mx.model <- omxAssignFirstParameters(mx.model)
  if (!is.null(RAM)) {
    mx.model <- .set.label.bounds(mx.model, RAM.lbound, RAM.ubound)
    ## Apply explicit overrides by label after OpenMx reconciles shared
    ## parameters, so an NA in another occurrence cannot erase the bound.
    mx.model <- .set.label.bounds(mx.model, A.lower.override, A.upper.override)
  }
  mx.model <- .tighten.feasibility(mx.model)
    
  warning.msg <- error.msg <- NULL
  if (run) {
    mx.fit <- tryCatch(
      withCallingHandlers(
        do.call(mxRun, c(list(mx.model, intervals=intervals,
                              suppressWarnings=suppressWarnings,
                              silent=silent), mxRun.Args)),
        warning=function(w) {
          warning.msg <<- c(warning.msg, conditionMessage(w))
          invokeRestart("muffleWarning")
        }),
      error=function(e) {
        error.msg <<- conditionMessage(e)
        e
      })
    if (inherits(mx.fit, "error")) stop(mx.fit)
    ## Retry a fit without status 0 or 1 (see .retry.bad.status()).
    status <- NULL
    if (!.status.ok(mx.fit)) {
      status <- mx.fit@output$status$code
      record <- function(expr) withCallingHandlers(
        tryCatch(expr, error=function(e) e),
        warning=function(w) {
          warning.msg <<- c(warning.msg, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      mx.fit <- .retry.bad.status(mx.fit, record=record,
                                  rerun=if (intervals) function(fit)
                                    do.call(mxRun, c(list(fit, intervals=TRUE,
                                                          suppressWarnings=TRUE,
                                                          silent=silent), mxRun.Args)))
    }
    ## Warnings are recorded so that they cannot abort the fit, and are
    ## shown after it unless suppressWarnings=TRUE.
    if (!suppressWarnings)
      for (w in unique(warning.msg)) warning(w, call.=FALSE)
    mx.fit <- .check.constraints(mx.fit, intervals=intervals)
    .warn.bad.status(mx.fit, status)
    if (!is.null(mx.fit$Tau2))
      .warn.Tau.not.psd(mxEvalByName("Tau2", mx.fit, compute=TRUE), "Tau2")
  } else {
    mx.fit <- mx.model
  }

  out <- list(call=match.call(), 
              Mmatrix=Mmatrix, 
              Tmatrix=Tmatrix,
              Vmatrix=Vmatrix, 
              Jmatrix=Jmatrix, 
              data=data,
              labels=list(obslabels=obslabels, ylabels=ylabels, vlabels=vlabels),
              mxModel.Args=mxModel.Args,
              mxRun.Args=mxRun.Args,
              subset.variables=subset.variables,
              subset.rows=subset.rows,
              #mx.model=mx.model,
              mx.fit=mx.fit,
              warning=warning.msg,
              error=error.msg)
  class(out) <- 'osmasem'
  out
}

coef.osmasem <- function(object, select=c("fixed", "all", "random"), ...) {
  if (!inherits(object, "osmasem"))
    stop("\"object\" must be an object of class \"osmasem\".")
    
  mx.coef <- coef(object$mx.fit)
  ## index for untransformed random effects (not the correct ones!) 
  index <- names(mx.coef) %in% .osmasem.random.labels(object)
  
  select <- match.arg(select)
  switch(select,
         fixed =  mx.coef <- mx.coef[!index],
         random = mx.coef <- mx.coef[index])
  if (select!="fixed") warning("\"Tau1_xx\" is not the variance component of the random effects.\nPlease use VarCorr() to get the variance component.\n")    
  mx.coef    
}

## Heterogeneity labels come from the actual random-effects matrices, not
## from substrings in arbitrary user-supplied fixed-effect labels.
.osmasem.random.labels <- function(object) {
  matrices <- object$Tmatrix
  if (is.null(matrices)) return(character())
  unique(unlist(lapply(matrices, function(m) {
    if (inherits(m, "MxMatrix")) m$labels[m$free] else NULL
  }), use.names=FALSE))
}

## Set lower and upper bounds (named lists keyed by parameter label) on the
## free parameters of an osmasem() model, in whichever matrix they are.
## Parameters that are not free in the model (e.g., error variances that
## osmasem() computes from the other parameters) are ignored with a warning.
.set.label.bounds <- function(model, lbound=NULL, ubound=NULL) {
  free <- names(omxGetParameters(model))
  ignored <- setdiff(names(c(lbound, ubound)), free)
  if (length(ignored))
    warning("Bounds on parameter(s) ", paste(unique(ignored), collapse=", "),
            " are ignored: they are not free parameters in the osmasem() model.",
            call.=FALSE)
  lbound <- lbound[intersect(names(lbound), free)]
  ubound <- ubound[intersect(names(ubound), free)]
  if (length(lbound))
    model <- omxSetParameters(model, labels=names(lbound), lbound=unlist(lbound))
  if (length(ubound))
    model <- omxSetParameters(model, labels=names(ubound), ubound=unlist(ubound))
  model
}

## Covariance of the estimates of osmasem() and osmasem2() fits. With
## mxConstraints, OpenMx's covariance (also the bread of imxRobustSE()) is
## replaced by .model.vcov() and the sandwich of .robust.observation.vcov().
.osmasem.vcov <- function(model, robust=FALSE) {
  if (!.has.constraints(model)) {
    if (robust) suppressMessages(imxRobustSE(model, details=TRUE))$cov
    else vcov(model)
  } else if (robust) {
    .robust.observation.vcov(model)
  } else {
    .model.vcov(model)
  }
}

vcov.osmasem <- function(object, select=c("fixed", "all", "random"),
                         robust=FALSE, ...) {
  if (!inherits(object, "osmasem"))
    stop("\"object\" must be an object of class \"osmasem\".")
  
  # labels of the parameters    
  ## my.name <- summary(object$mx.fit)$parameters$name
  my.name <- names(omxGetParameters(object$mx.fit))
  my.name <- my.name[!is.na(my.name)]

  ## index for untransformed random effects (not the correct ones!) 
  index.random <- my.name %in% .osmasem.random.labels(object)

  select <- match.arg(select)
  switch(select,
    fixed = my.name <- my.name[!index.random],
    random = my.name <- my.name[index.random])

  out <- .osmasem.vcov(object$mx.fit, robust=robust)
  if (select!="fixed") warning("\"Tau1_xx\" is not the variance component of the random effects.")
  out[my.name, my.name, drop=FALSE]
}



#' Extract Variance-Covariance Matrix of the Random Effects
#' 
#' It extracts the variance-covariance matrix of the random effects (variance
#' component) from either the \code{meta} or \code{osmasem} objects.
#' 
#' 
#' @param x An object returned from either class \code{meta} or \code{osmasem}
#' @param \dots Further arguments; currently none is used
#' @return A variance-covariance matrix of the random effects.
#' @note It is similar to \code{coef(object, select="random")} in tssem. The
#' main difference is that \code{coef()} returns a vector while
#' \code{VarCorr()} returns its correspondent matrix.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{coef}}, \code{\link[metaSEM]{vcov}}
#' @keywords methods
#' @examples
#' 
#' ## Multivariate meta-analysis on the log of the odds
#' ## The conditional sampling covariance is 0
#' bcg <- meta(y=cbind(ln_Odd_V, ln_Odd_NV), data=BCG,
#'             v=cbind(v_ln_Odd_V, cov_V_NV, v_ln_Odd_NV))
#' VarCorr(bcg)
#' 
VarCorr <- function(x, ...) {
  if (all(!is.element(c("meta", "osmasem", "osmasem3L"), class(x))))
    stop("\"x\" must be an object of either class \"meta\", \"osmasem\" or \"osmasem3L\".")

  ## switch(class(x)[1],
  ## Assuming only one class in x
  switch(class(x)[class(x) %in% c("meta", "osmasem", "osmasem3L")],
         meta = out <- mxEvalByName("Tau", x$mx.fit),
         osmasem = {
           out <- mxEvalByName("Tau2", x$mx.fit)
           
           ## ## No. of variables in the model
           ## p <- ncol(out)
           ## names.tau2 <- paste0("Tau2_", seq_len(p))

           ## Get the names Tau1_1, Tau1_2... or Tau1_a, Tau1_b... from x$Tmatrix
           ## Replace Tau1 with Tau2
           names.tau2 <- gsub("Tau1", "Tau2", c(x$Tmatrix$vecTau1$labels))
           dimnames(out) <- list(names.tau2, names.tau2)},
         osmasem3L = {
           ## Tau2 (within)
           Tau2W <- mxEvalByName("Tau2W", x$mx.fit)
           ## Labels ended with "W"
           names.tau2W <- gsub("W$", "", c(x$TmatrixW$vecTau1$labels))
           names.tau2W <- gsub("Tau1", "Tau2", names.tau2W)               
           dimnames(Tau2W) <- list(names.tau2W, names.tau2W)

           ## Tau2 (between)
           Tau2B <- mxEvalByName("B.Tau2B", x$mx.fit)
           ## Labels ended with "B" 
           names.tau2B <- gsub("B$", "", c(x$TmatrixB$vecTau1$labels))
           names.tau2B <- gsub("Tau1", "Tau2", names.tau2B)
           dimnames(Tau2B) <- list(names.tau2B, names.tau2B)
           out <- list(Tau2W=Tau2W, Tau2B=Tau2B)},
         stop("Invalid class:", class(x)))

  out    
}


## Fit a saturated model with either a diagonal or symmetric variance component of random effects
.osmasemSatIndMod <- function(osmasem.obj=NULL, model=c("Saturated", "Independence"),
                              Std.Error=FALSE, Tmatrix, data, 
                              mxModel.Args=NULL, mxRun.Args=NULL, suppressWarnings=TRUE,
                              silent=FALSE, run=TRUE, ...) {

  if (!is.null(osmasem.obj)) {
    if (!is.element("osmasem", class(osmasem.obj)))
      stop("\"osmasem.obj\" must be an object of class \"osmasem\".")
    
    if (missing(Tmatrix)) {
      Tmatrix <- osmasem.obj$Tmatrix
    }                
    ## No. of variables in the model
    p <- nrow(Tmatrix$Cor$values)
    data <- osmasem.obj$data
    ylabels <- osmasem.obj$labels$ylabels
    Vmatrix <- osmasem.obj$Vmatrix        
    subset.rows <- osmasem.obj$subset.rows
  } else {
    ## Build a saturated/independence model directly from Cor2DataFrame data.
    if (missing(data) || is.null(data$data) || is.null(data$ylabels) ||
        is.null(data$vlabels))
      stop("'data' must be an object returned by Cor2DataFrame().")
    if (missing(Tmatrix))
      stop("'Tmatrix' is required when 'osmasem.obj' is NULL.")
    ylabels <- data$ylabels
    p <- length(ylabels)
    ## Create known sampling variance covariance matrix
    Vmatrix <- create.V(data$vlabels, type="Symm", as.mxMatrix=TRUE)
    subset.rows <- rep(TRUE, nrow(data$data))
  }    

  model <- match.arg(model)
  
  if (model=="Saturated") {
    ## Saturated model
    st.values <- apply(data$data[, ylabels, drop=FALSE], 2,  mean, na.rm=TRUE)
    Mmatrix <- mxMatrix(type="Full", free=TRUE, values=st.values,
                        labels=paste0("Mu", seq_len(p)),
                        nrow=1, ncol=p, dimnames=list(NULL, ylabels), name="Mu")
    
  } else {
    ## Independence model
    Mmatrix <- mxMatrix(type="Full", free=FALSE, labels=paste0("Mu", seq_len(p)),
                        nrow=1, ncol=p, name="Mu")        
  }

  ## Filter the data for the analysis
  mx.data <- data$data[subset.rows, , drop=FALSE]
  
  ## Dirty trick to avoid the warning of "no visible binding for global variable"
  mx.model <- eval(parse(text="mxModel(model=model,
                               mxData(observed=mx.data, type='raw'),
                               mxExpectationNormal(covariance='expCov',
                                                   means='Mu',
                                                   dimnames=ylabels),
                               mxFitFunctionML(), Mmatrix, Tmatrix, Vmatrix,
                               mxAlgebra(Tau2+V, name='expCov'))"))
    
  ## Add additional arguments to mxModel
  if (!is.null(mxModel.Args)) {
    for (i in seq_along(mxModel.Args)) {
      mx.model <- mxModel(mx.model, mxModel.Args[[i]])
    }
  }

  ## It may speed up the analysis
  if (Std.Error==FALSE) {
    mx.model <- mxOption(mx.model, "Calculate Hessian", "No")
    mx.model <- mxOption(mx.model, "Standard Errors", "No")
  }
  
  ## Return mx model without running the analysis
  if (run==FALSE) {
    return(mx.model)
  } else {
    mx.fit <- tryCatch(do.call(mxRun, c(list(mx.model, silent=silent,
                                             suppressWarnings=suppressWarnings),
                                        mxRun.Args)), error = function(e) e)
  }
    
  # try to run it with error message as output
  ## if (inherits(mx.fit, "error")) {
  ##     warning(print(mx.fit))
  ## }
  mx.fit
}

anova.osmasem <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$mx.fit)
  comparison <- lapply(list(...), function(x) x$mx.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

summary.osmasem <- function(object, fitIndices=FALSE, numObs,
                            robust=FALSE, ...) {
  if (!inherits(object, "osmasem"))
    stop("\"object\" must be an object of class \"osmasem\".")
  
  ## Ind.stat <- summary(.osmasemSatIndMod(object, model="Independence",
  ##                                    Std.Error=FALSE, silent=TRUE))   

  ## If numObs is not provided, use the total N
  if (missing(numObs)) {
    ## numObs <- sum(object$data$n)-length(object$data$n)
    numObs <- sum(object$data$n[object$subset.rows])
  }
  
  ## Calculate chi-square statistic and other fit indices
  if (fitIndices) {   ## Check if the dimensions of Tmatrix and Vmatrix match.
    ## They can be different for rcmasem
    ## No. of effect sizes
    p <- dim(object$Vmatrix$values)[1]
    
    if (length(object$Tmatrix$vecTau1$labels) != p) {   ## Check if there are correlations,
      ## If yes -> Symm, no -> Diag
      ## if (sum(object$Tmatrix$Cor$free)>0) RE.type="Symm" else RE.type="Diag"
      ## T0 <- create.Tau2(no.var=p, RE.type=RE.type)

      ## TODO: For rcmasem, the variance component cannot be Diag;
      ## otherwise, the saturated model fits poorer than the rcmasem.
      T0 <- create.Tau2(no.var=p, RE.type="Symm")
      
      Sat.stat <- .osmasemSatIndMod(object, model="Saturated", Tmatrix=T0,
                                    Std.Error=FALSE, silent=TRUE)
      Ind.stat <- .osmasemSatIndMod(object, model="Independence", Tmatrix=T0,
                                    Std.Error=FALSE, silent=TRUE)            
    } else {
      Sat.stat <- .osmasemSatIndMod(object, model="Saturated",
                                    Std.Error=FALSE, silent=TRUE)
      Ind.stat <- .osmasemSatIndMod(object, model="Independence", 
                                    Std.Error=FALSE, silent=TRUE) 
    }
       
    ## Sat.model=FALSE for saturated model if there are either errors or nonconvergent.
    if ( inherits(Sat.stat, "error") | !(Sat.stat$output$status$code %in% c(0,1)) ) {
      Sat.model <- FALSE
    } else {
      Sat.stat <- summary(Sat.stat)
      if (Sat.stat$degreesOfFreedom <= 0) {
        Sat.model <- FALSE
      } else {
        Sat.model <- TRUE
      }
    }

    ## Ind.model=FALSE for independence model if there are either errors or nonconvergent.
    ## NA is acceptable for independence model as it means that optimization was not attempted.
    if ( inherits(Ind.stat, "error") | !(Ind.stat$output$status$code %in% c(NA, 0, 1)) ) {
      Ind.model <- FALSE
    } else {
      Ind.stat <- summary(Ind.stat)
      if (Ind.stat$degreesOfFreedom <= 0) { 
        Ind.model <- FALSE
      } else {
        Ind.model <- TRUE
      }
    }

    ## If sat.model is defined
    if (Sat.model) {
      if (Ind.model) {
        out <- summary(object$mx.fit,
                       SaturatedLikelihood=Sat.stat$Minus2LogLikelihood,
                       SaturatedDoF=Sat.stat$degreesOfFreedom,
                       IndependenceLikelihood=Ind.stat$Minus2LogLikelihood,
                       IndependenceDoF=Ind.stat$degreesOfFreedom,
                       numObs=numObs, ...)
      } else {
        out <- summary(object$mx.fit,
                       SaturatedLikelihood=Sat.stat$Minus2LogLikelihood,
                       SaturatedDoF=Sat.stat$degreesOfFreedom,
                       numObs=numObs, ...)
        warning("There are errors in either fitting the independence model or its degree of freedom being non-positive.\n")
      }
      ## if Sat.model is not defined, no chi-square and fit indices    
    } else {
      out <- summary(object$mx.fit, numObs=numObs, ...)
      warning("There are errors in either fitting the saturated model or its degree of freedom being non-positive.\n")
    }
    ## fitIndices=FALSE 
  } else {
    out <- summary(object$mx.fit, numObs=numObs, ...)
  }

  ## Modified the SE with robust SE, or with model-based SE when there are
  ## mxConstraints
  if (robust || .has.constraints(object$mx.fit)) {
    SE <- sqrt(diag(.osmasem.vcov(object$mx.fit, robust=robust)))
    out$parameters[, "Std.Error"] <- SE[out$parameters$name]
  }

  ## Additional output to the mx summary
  out$parameters$`z value` <- out$parameters$Estimate / out$parameters$Std.Error
  out$parameters$`Pr(>|z|)` <- 2*pnorm(-abs(out$parameters$`z value`))
  out
}



#' Calculate the R2 in OSMASEM and OSMASEM3L
#' 
#' It calculates the R2 of the moderators in explaining the variances in the
#' heterogeneity variances.
#' 
#' 
#' @param model1 An object inheriting from class \code{osmasem}.
#' @param model0 An object inheriting from class \code{osmasem}.
#' @param R2.truncate Whether to truncate the negative R2 to zero.
#' @return \code{model1} and \code{model0} are the models with and without the
#' moderators, respectively. The function does not check whether the models are
#' nested. It is the users' responsibility to make sure that the models with
#' and without the moderators are nested. It returns a list of the diagonals of
#' the heterogeneity variances of the models without and with the moderators,
#' and the R2.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{osmasem}}
#' @keywords osmasem osmasem3L
osmasemR2 <- function(model1, model0, R2.truncate=TRUE) {
  if (!inherits(model0, "osmasem") || !inherits(model1, "osmasem"))
    stop("Both \"model0\" and \"model1\" must be objects of class \"osmasem\".")

  ## Tau2.0 <- diag(eval(parse(text = "mxEval(Tau2, model0$mx.fit)")))       
  ## Tau2.1 <- diag(eval(parse(text = "mxEval(Tau2, model1$mx.fit)")))
  Tau2.0 <- diag(VarCorr(model0))
  Tau2.1 <- diag(VarCorr(model1))               

  R2 <- (Tau2.0-Tau2.1)/Tau2.0

  ## Truncate negative R2 to 0
  if (R2.truncate) R2[R2<0] <- 0
    
  p <- nrow(model0$Tmatrix$Cor$values)        
  names.tau2 <- paste0("Tau2_", seq_len(p), "_", seq_len(p))
  names(R2) <- names(Tau2.0) <- names(Tau2.1) <- names.tau2
    
  list(Tau2.0=Tau2.0, Tau2.1=Tau2.1, R2=R2)
}



#' Calculate the SRMR in OSMASEM and OSMASEM3L
#' 
#' It calculates the standardized root mean squared residuals (SRMR) in OSMASEM
#' and OSMASEM3L.
#' 
#' 
#' @param x An OSMASEM object without any moderators.
#' @return It calculates the model implied correlation matrix and its saturated
#' counterpart to calculate the SRMR. It should be noted that the heterogeneity
#' variances are ignored in the calculations.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{osmasem}}, \code{\link[metaSEM]{Nohe15}}
#' @keywords osmasem osmasem3L
osmasemSRMR <- function(x) {
  ## Check if there are moderators in either A1 or S1
  if (!is.null(x$Mmatrix$A1) | !is.null(x$Mmatrix$S1))
    stop("Moderators are not allowed in calculating the SRMR in OSMASEM.\n")
  
  ## Saturated model which should be very close to the stage 1 results in the TSSEM
  Sat.stat <- .osmasemSatIndMod(x, model="Saturated", Std.Error=FALSE, silent=TRUE)

  ## Similar to the sample correlation matrix
  sampleR <- vec2symMat(mxEvalByName("Mu", Sat.stat), diag=FALSE)

  ## Model implied correlation matrix
  impliedR <- mxEvalByName("impliedR", x$mx.fit)

  ## SRMR (Cheung, 2015 book, Eq. 2.36)
  sqrt(mean(vechs(sampleR-impliedR)^2))
}

## RAMmodelV <- function(Amatrix, Smatrix, Fmatrix, Mmatrix, Vmatrix,
##                       data, obs.var, ...) {

##     Sfull <- mxAlgebra(Smatrix + Vmatrix, name="Sfull")
##     fitfun <- mxExpectationRAM(A="Amatrix", S="Sfull", F="Fmatrix",
##                                M="Mmatrix", dimnames=obs.var)
##     fit <- mxFitFunctionML()
##     mxModel(mxData(observed=data, type="raw"),
##             Amatrix, Smatrix, Fmatrix, Mmatrix, Vmatrix, Sfull,
##             fitfun, fit, ...)
## }

## NORMmodelV <- function(Mu, Sigma, Fmatrix, Vmatrix, data, obs.var, ...) {

##     Sfull <- mxAlgebra(Sigma + Vmatrix, name="Sfull")
##     fitfun <- mxExpectationNormal(covariance="Sfull",
##                                   means=Mu@name,
##                                   dimnames=obs.var)
##     fit <- mxFitFunctionML()
##     mxModel(mxData(observed=data, type="raw"),
##             Mu, Sigma, Fmatrix, Vmatrix, Sfull, fitfun, fit, ...)
## }
