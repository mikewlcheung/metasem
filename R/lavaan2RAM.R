#' Convert \code{lavaan} models to RAM models
#' 
#' It converts models specified in \code{lavaan} model syntax to RAM models.
#' 
#' It uses the \code{\link[lavaan]{model.syntax}} to do the conversion.
#' Experimental: functions of parameters (`:=` in lavaan) and constraints
#' (`==`, `>`, and `<` in lavaan) will be converted to mxAlgebra and
#' mxConstraint in OpenMx. As there are differences between lavaan and OpenMx,
#' they may not work properly.
#' 
#' @param model A character string of model using the lavaan model syntax.
#' @param obs.variables A character vector of the observed variables. The
#' observed variables in the RAM specification will follow the order specified
#' in \code{obs.variables}. It is important to check whether the order of the
#' observed variables matches the order in the dataset.
#' @param A.notation A character string to be used in the A matrix if the
#' labels are not included in the lavaan model. For example, the label will be
#' "yONx" for regressing "y" on "x".
#' @param S.notation A character string to be used in the S matrix if the
#' labels are not included in the lavaan model. For example, the label will be
#' "yWITHx" for the covariance between "y" with "x" and "yWITHy" for the
#' (error) variance of "y".
#' @param M.notation A character string to be used in the M matrix if the
#' labels are not included in the lavaan model. For example, the label will be
#' "ymean" for the mean of "y" if \code{M.notation="mean"}.
#' @param A.start A numeric value of starting value for the \code{Amatrix} when
#' the starting values are not provided.
#' @param S.start A numeric value of starting value for the \code{Smatrix} when
#' the starting values are not provided.
#' @param M.start A numeric value of starting value for the \code{Mmatrix} when
#' the starting values are not provided.
#' @param auto.var Logical. If \code{TRUE}, the residual variances and the
#' variances of exogenous latent variables are included in the model and set
#' free. See \code{\link[lavaan]{model.syntax}}.
#' @param std.lv Logical. If \code{TRUE} (the default), the metric of each
#' latent variable is determined by fixing their variances to 1.0. If
#' \code{FALSE}, the metric of each latent variable is determined by fixing
#' the factor loading of the first indicator to 1.0. See
#' \code{\link[lavaan]{model.syntax}}. \code{std.lv=TRUE} always fixes every
#' latent variable's variance at 1, even overriding an explicit \code{"~~"}
#' specification for it in \code{model} (e.g. \code{"f ~~ c(vf1,vf2)*f"} to
#' free a factor's variance per group) with no error -- the label is simply
#' discarded, since the parameter is fixed rather than estimated. Pass
#' \code{std.lv=FALSE} whenever \code{model} identifies a latent variable's
#' scale by fixing one of its loadings instead (e.g. \code{"f =~ 1*y1"} or
#' \code{"f =~ c(1,1)*y1"} for a multiple-group model), which is
#' incompatible with also fixing its variance. \code{lavaan2RAM()} warns
#' when this happens (any parameter given an explicit label in \code{model}
#' that ends up fixed rather than free), but the warning is general --
#' checking \code{std.lv} first is faster than working backward from it.
#' @param ngroups Number. The number of groups in the \code{model}. See
#' \code{\link[lavaan]{model.syntax}}. Parameters without an explicit label
#' in \code{model} are free and distinct in each group by default (as in
#' \code{lavaan}); use \code{lavaan}'s \code{c(label, label)} syntax to
#' constrain a parameter to be equal across groups. Has no effect (and is
#' rejected with an error) on two-level (\code{level:}) models -- see below.
#' @param \dots Further arguments to be passed to
#' \code{\link[lavaan]{model.syntax}}
#' @return A list of RAM specification with \code{A}, \code{S}, \code{F}, and
#' \code{M} matrices. When \code{ngroups > 1}, a list of such lists (one per
#' group) is returned instead, with class \code{"RAM_multigroup"}. When
#' \code{model} uses \code{lavaan}'s two-level syntax (\code{level: 1} /
#' \code{level: 2} blocks), a list with \code{within} and \code{between}
#' elements (each themselves a list of \code{A}, \code{S}, \code{F}, and
#' \code{M} matrices) is returned instead, with class \code{"RAM_twolevel"}.
#' Level-1 (within) manifest variables default to a fixed-at-zero mean (the
#' overall intercept lives at the between level by convention); level-2
#' (between) "shadow" variables backing a level-1 manifest variable default
#' to a free mean, while genuine level-2-only latent factors default to
#' fixed-at-zero, mirroring the usual latent-variable convention. Genuine
#' level-2-only observed covariates (e.g. a real cluster-level covariate
#' with its own data column, as opposed to a level-1 manifest variable's
#' between-cluster component) are not yet supported and raise an error,
#' rather than being silently treated as an ordinary latent variable with
#' no connection to any data; neither is combining two-level with multiple
#' groups (also an error).
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[semPlot]{ramModel}}, \code{\link[metaSEM]{Becker92}},
#' \code{\link[metaSEM]{Becker09}}, \code{\link[metaSEM]{Digman97}},
#' \code{\link[metaSEM]{Hunter83}}, \code{\link[metaSEM]{as.mxMatrix}},
#' \code{\link[metaSEM]{checkRAM}}
#' @keywords methods tssem
#' @examples
#' 
#' ## Regression model on correlation matrix
#' model1 <- "## y is modelled by x1, x2, and x3
#'            y ~ b1*x1 + b2*x2 + b3*x3
#'            ## Fix the independent variables at 1
#'            x1 ~~ 1*x1
#'            x2 ~~ 1*x2
#'            x3 ~~ 1*x3
#'            ## Declare the correlations among the independent variables
#'            x1 ~~ x2
#'            x1 ~~ x3
#'            x2 ~~ x3
#'            ## Constraint
#'            b3 == b1 + b2
#'            ## Function of parameters
#'            fn1 := b1*b2^b3"
#' 
#' ## Compare the arrangements of variables with and without
#' ## specifying the obs.variables arguments. 
#' lavaan2RAM(model1, obs.variables=c("y", "x1", "x2", "x3"))
#' 
#' ## Two-factor CFA model
#' model2 <- "f1 =~ x1 + x2 + x3
#'            f2 =~ x4 + x5 + x6
#'            ## Declare the correlation between f1 and f2
#'            ## and label it with cor_f1f2
#'            f1 ~~ cor_f1f2*f2"
#' 
#' lavaan2RAM(model2)
#' 
#' ## Regression model with the mean structure
#' model3 <- "y ~ x
#'            ## Intercept of y
#'            y ~ 1
#'            ## Mean of x
#'            x ~ 1"
#' 
#' lavaan2RAM(model3)
#' 
#' ## Multiple group regression model
#' ## Different intercepts with a common slope
#' model4 <- "y ~ c(a1, a2)*1 + c(b, b)*x"
#' 
#' lavaan2RAM(model4, ngroups=2)
#'
#' ## Two-level random-intercept regression model
#' model5 <- "level: 1
#'              y ~ b1*x
#'              y ~~ residW*y
#'            level: 2
#'              y ~ 1
#'              y ~~ residB*y"
#'
#' lavaan2RAM(model5)
#'
lavaan2RAM <- function(model, obs.variables = NULL, A.notation="ON", S.notation="WITH",
                       M.notation="mean", A.start=0.1, S.start=0.5, M.start=0,
                       auto.var = TRUE, std.lv = TRUE,
                       ngroups=1, ...) {
  ## if (!requireNamespace("lavaan", quietly=TRUE))    
  ##     stop("\"lavaan\" package is required for this function.")

  ## Default: fix the latent independent variables at 1
  my.model <- lavaan::lavaanify(model, fixed.x = FALSE, auto.var=auto.var,
                                std.lv=std.lv, ngroups=ngroups, ...)

  ## Warn about a parameter that the user gave an explicit label but that
  ## lavaanify() has fixed anyway (free==0): a label on a fixed value does
  ## nothing (there is no free parameter for it to equate or report), so
  ## this is virtually always unintentional. By far the most common cause
  ## is a latent variable's own "~~" variance/covariance clashing with
  ## std.lv=TRUE (this function's default), which always fixes every latent
  ## variable's variance at 1 for identification, silently discarding any
  ## explicit "~~" specification for it with no error -- e.g. writing
  ## "f ~~ c(vf1,vf2)*f" to free a factor's variance per group has no effect
  ## unless std.lv=FALSE is also passed (needed anyway here, since fixing a
  ## loading, e.g. "f =~ 1*y1", for identification is incompatible with also
  ## fixing the variance via std.lv=TRUE). Not specific to std.lv, though --
  ## the same "labelled but fixed" signature also flags e.g. a covariate
  ## mistakenly labelled under fixed.x=TRUE elsewhere. label=="" (not NA)
  ## for an unlabelled row, confirmed against lavaanify()'s own output.
  ## Excludes ":="/"=="/"<"/">" rows -- these are lavaan's own defined-
  ## parameter/constraint operators, not ordinary model parameters: a ":="
  ## row is *always* free==0 by design (it is a deterministic mxAlgebra of
  ## other parameters, not something ML estimates directly) and its
  ## "label" column holds the defined parameter's own name (e.g. "fn1" for
  ## "fn1 := b1+b2"), not a constraint label -- confirmed against
  ## lavaanify()'s own output, not assumed; every ":=" row would otherwise
  ## false-positive here.
  fixed.but.labeled <- my.model$free==0 & nzchar(my.model$label) &
    !(my.model$op %in% c(":=", "==", "<", ">"))
  if (any(fixed.but.labeled)) {
    bad <- my.model[fixed.but.labeled, ]
    msg <- paste0(bad$lhs, " ", bad$op, " ", bad$rhs, " (label \"", bad$label, "\")")
    warning("The following parameter(s) were given an explicit label in ",
           "'model' but ended up FIXED (not free) after lavaanify()'s own ",
           "defaults were applied, so the label has no effect: ",
           paste(msg, collapse="; "), ". This most commonly happens for a ",
           "latent variable's own variance/covariance under std.lv=TRUE ",
           "(the default here), which always fixes it at 1 regardless of ",
           "any explicit \"~~\" specification -- pass std.lv=FALSE if you ",
           "intended it to be free (e.g. because the model already fixes a ",
           "loading for identification instead).", call.=FALSE)
  }

  ## Two-level ("level: 1" / "level: 2") syntax: lavaanify() reports this
  ## via a "level" column instead of the usual "group" column (the two are
  ## mutually exclusive unless nested "group:"/"level:" block syntax is
  ## used, which is not yet supported here -- see dev/PLAN-two-level-
  ## multigroup-sem.md, phase 3).
  if ("level" %in% colnames(my.model)) {
    if ("group" %in% colnames(my.model)) {
      stop("Combined two-level and multiple-group models (nested \"group:\"",
          "/\"level:\" syntax) are not yet supported by lavaan2RAM().\n")
    }
    if (ngroups > 1) {
      stop("\"ngroups\" > 1 has no effect on two-level (\"level:\") models; ",
          "combined two-level and multiple-group models are not yet ",
          "supported by lavaan2RAM().\n")
    }
    return(.lavaan2RAM.twolevel(my.model, obs.variables=obs.variables,
                                A.notation=A.notation, S.notation=S.notation,
                                M.notation=M.notation, A.start=A.start,
                                S.start=S.start, M.start=M.start))
  }

  ## Maximum no. of groups
  max.gp <- max(my.model$group)

  ## Empty list to store the matrices per group
  out <- list()
    
  for (gp in seq_len(max.gp)) {
    ## Select the ith group
    mod <- my.model[my.model$group==gp, ]

    ## if (any...) is required to avoid error when there is no element for the assignment
    ## set the starting values in A if NA
    if (any((mod$op=="=~"|mod$op=="~")&is.na(mod$ustart))) {            
      mod[(mod$op=="=~"|mod$op=="~")&is.na(mod$ustart), ]$ustart <- A.start
    }
        
    ## set the starting values in M if NA
    if (any(mod$op=="~1"&is.na(mod$ustart))) {
      mod[mod$op=="~1"&is.na(mod$ustart), ]$ustart <- M.start      
    }
        
    ## set the starting values in S and free parameters if NA (variances)
    if (any(mod$op=="~~"&is.na(mod$ustart)&(mod$lhs==mod$rhs))) {
      mod[mod$op=="~~"&is.na(mod$ustart)&(mod$lhs==mod$rhs), ]$ustart <- S.start
    }
        
    ## Set the starting values in S and free parameters if NA (covariances)
    if (any(mod$op=="~~"&is.na(mod$ustart)&(mod$lhs!=mod$rhs))) {
      mod[mod$op=="~~"&is.na(mod$ustart)&(mod$lhs!=mod$rhs), ]$ustart <- 0
    }

    ## all variables
    ## Removed sort(); otherwise, the variables will be arranged as x1, x10, x2, x3...
    all.var <- unique(c(mod$lhs, mod$rhs))
        
    ## latent variables: (with indicators)
    latent <- unique(mod[mod$op== "=~", ]$lhs)
    ## observed variables: not latent
    observed <- all.var[!(all.var %in% latent)]
    ## remove empty string "" when there are mean structure
    observed <- observed[observed !=""]
        
    ## check whether observed in model = observed in argument
    if (!is.null(obs.variables)) {
      if (!identical(sort(observed), sort(obs.variables))) {
        stop("Names in \"obs.variables\" do not agree with those in model.\n")
      } else {
        ## arrange the observed variables according to obs.var argument
        observed <- obs.variables
      }
    }

    ## if there are latent variables
    if (length(latent)>0) {
      ## arrange variable list as observed + latent
      all.var <- c(observed, latent)
    } else {
      all.var <- observed
    }
    
    no.lat <- length(latent)
    no.obs <- length(observed)
    no.all <- no.lat+no.obs
    
    Amatrix <- matrix(0, ncol=no.all, nrow=no.all, dimnames=list(all.var, all.var))
    Smatrix <- matrix(0, ncol=no.all, nrow=no.all, dimnames=list(all.var, all.var))
    ## Fixed Mmatrix by setting the default to estimate the means;
    ## otherwise, the mean structure is misspecified.
    ## Defaults: observed variables are free, whereas latent variables are fixed at 0
    ## This default label bypasses the per-row labelling loop below (there is
    ## no explicit "~1" row to attach a group suffix to), so the group
    ## suffix must be added here directly for the same reason as in that
    ## loop: keep auto-generated labels distinct per group by default.
    mean.suffix <- if (max.gp > 1) paste0("_g", gp) else ""
    Mmatrix <- matrix(c(paste0("0*", observed, "mean", mean.suffix), rep("0", no.lat)),
                      nrow=1, ncol=no.all, dimnames=list(1, all.var))

    ## Prepare the labels
    for (i in seq_len(nrow(mod))) {
      ## if there is no label
      if (mod[i, ]$label=="") {
        switch(mod[i, ]$op,
               "=~" = mod[i, ]$label <- paste0(mod[i, ]$rhs, A.notation, mod[i, ]$lhs),
               "<~" = mod[i, ]$label <- paste0(mod[i, ]$lhs, A.notation, mod[i, ]$rhs),
               "~"  = mod[i, ]$label <- paste0(mod[i, ]$lhs, A.notation, mod[i, ]$rhs),
               "~~" = mod[i, ]$label <- paste0(mod[i, ]$lhs, S.notation, mod[i, ]$rhs),
               "~1" = mod[i, ]$label <- paste0(mod[i, ]$lhs, M.notation))
        ## Auto-generated labels must stay distinct per group by default,
        ## as lavaan itself does not equate unlabelled parameters across
        ## groups. Only labels the user explicitly shares across groups
        ## (e.g., c(b,b)*x) are equal; those already carry the same string
        ## in "label" before reaching this branch, so they are untouched.
        if (max.gp > 1) {
          mod[i, ]$label <- paste0(mod[i, ]$label, "_g", gp)
        }
      }
    }

    ## Replace NA to 0 in ustart if there are still NA
    ## mod$ustart[is.na(mod$ustart)] <- 0
        
    ## keys in as.mxMatrix format
    key <- with(mod, ifelse(free==0, yes=ustart, no=paste(ustart, label, sep="*")))  
        
    for (i in seq_len(nrow(mod))) {
      my.line <- mod[i, ]
      ## lavaanify() auto-adds a FIXED ("user"==0, "free"==0, "ustart"==0)
      ## "~1" row for any OBSERVED variable that never gets its own
      ## explicit "~1" line, whenever a meanstructure is implied elsewhere
      ## in the model (e.g. a DV's own explicit "~1", as in "y ~ a1*1").
      ## Left alone, this auto row overwrites the free-mean default set
      ## above (Mmatrix's "0*<name>mean" initialisation) with a hard,
      ## non-free 0 -- silently contradicting this function's own
      ## documented default ("observed variables are free"). Confirmed
      ## empirically, not just suspected: for a covariate with a genuinely
      ## nonzero population mean (e.g. an uncentered predictor), this
      ## corrupts more than just the mean -- the matching "~~" (variance)
      ## parameter is then fit against data implicitly assumed mean-zero,
      ## so it absorbs the UNCENTERED second moment (mean^2 + variance)
      ## instead of the actual variance, and -2LL/fit statistics diverge
      ## sharply from an equivalent direct lavaan fit. An EXPLICIT
      ## user-written "~1" line (fixing a mean at a specific value, or
      ## freeing it with its own label) has "user"==1 and is unaffected --
      ## only lavaanify's own still-at-its-default auto rows are skipped.
      ##
      ## EXCEPT for an observed variable that is itself a "=~" INDICATOR
      ## of some latent factor (e.g. "f =~ 1*yi", the single-indicator
      ## "marker" trick used throughout this package's own meta-analysis
      ## functions, going back to meta()/meta3L()/metaFIML()): freeing
      ## such an indicator's own mean on top of its latent factor's own
      ## free mean (e.g. "f ~ mu*1") makes the two literally redundant --
      ## only their SUM is identified, not each individually. Confirmed
      ## this is not just a metaSEM quirk: real lavaan::sem() has the
      ## EXACT same non-identification for this pattern (explicit
      ## "Could not compute standard errors ... not identified" warning
      ## on this package's own canonical "f =~ 1*yi; f ~ mu*1" example),
      ## so blindly matching lavaan's literal free-by-default behaviour
      ## here would just import lavaan's own edge-case problem into the
      ## single idiom this entire package is built around. Keeping such
      ## an indicator's own mean fixed at 0 by default (as before) sides
      ## with metaSEM's own established, pervasive, well-identified
      ## convention over literal lavaan parity in this one case.
      if (my.line$op=="~1" && my.line$user==0 && my.line$free==0 &&
          my.line$lhs %in% observed &&
          !(my.line$lhs %in% mod[mod$op=="=~", "rhs"])) next
      switch(my.line$op,
             ## lhs: IV; rhs: DV
             "=~" = Amatrix[my.line$rhs, my.line$lhs] <- key[i],
             ## lhs: DV; rhs: IV
             "<~" = Amatrix[my.line$lhs, my.line$rhs] <- key[i],
             "~"  = Amatrix[my.line$lhs, my.line$rhs] <- key[i],
             "~~" = Smatrix[my.line$lhs, my.line$rhs] <-
               Smatrix[my.line$rhs, my.line$lhs] <- key[i],
             ## means
             "~1" = Mmatrix[1, my.line$lhs] <- key[i]
             )  ## from switch
    }  ## from for loop

    Fmatrix <- create.Fmatrix(c(rep(1, no.obs), rep(0, no.lat)), as.mxMatrix=FALSE)
    dimnames(Fmatrix) <- list(observed, all.var)
    
    out[[gp]] <- list(A=Amatrix, S=Smatrix, F=Fmatrix, M=Mmatrix)
  }

  ## Add group names, 1, 2, 3... to the list
  names(out) <- seq_along(out)
    
  ## If there are constraints such as .p1.==.p2.; remove them first
  ## otherwise, .p1.==.p2. will create an empty list in mxalgebra
  if (length(grep("^\\.", my.model$lhs)) >0 ) {
    my.model <- my.model[-grep("^\\.", my.model$lhs), ]
  }
    
  ## Check if there are constraints or algebras
  if (any(my.model$group==0)) {
        
    ## An empty list to store mxAlgebra and mxConstraint
    mxalgebra <- list()
    con_index <- 1   # A counter for no. of constraints, started from 1

    ## Constraints and algebras only
    y <- my.model[my.model$group==0, , drop=FALSE]
    
    for (i in seq_len(nrow(y))) {
      switch(y[i, 'op'],
             ## mxAlgebra
             ":=" = { eval(parse(text=paste0(y[i,'lhs'],
                                             "<- mxAlgebra(", y[i,'rhs'],
                                             ", name=\"", y[i,'lhs'], "\")")))
                                             eval(parse(text=paste0("mxalgebra <-
                                                  c(mxalgebra, ", y[i, 'lhs'],
                                                  "=", y[i, 'lhs'], ")")))  },
             ## Default condition to test if there are mxConstraints
             if (y[i, 'op'] %in% c("==", ">", "<")) {
               eval(parse(text=paste0("constraint", con_index, " <- mxConstraint(", y[i, 'lhs'],
                                      y[i, 'op'], y[i, 'rhs'], ", name=\"constraint", con_index, "\")")))
               eval(parse(text=paste0("mxalgebra <- c(mxalgebra, constraint",
                                      con_index, "=constraint", con_index, ")")))
               con_index <- con_index + 1
             }
             )  ## from switch
    }  ## from for loop

    ## Append mxalgebra to out[[1]]
    out[[1]]  <- list(A=out[[1]]$A, S=out[[1]]$S, F=out[[1]]$F, M=out[[1]]$M,
                      mxalgebras=mxalgebra)
  }  ## End if
  ## else {
  ##     out[[1]]  <- list(A=out[[1]]$A, S=out[[1]]$S, F=out[[1]]$F, M=out[[1]]$M)
  ## }


  ## Output the first list instead of a list of one item when there is only 1 group
  if (max.gp==1) {
    out  <- out[[1]]
  } else {
    ## Tag multi-group output so that sem() can dispatch on class() rather
    ## than sniffing the list structure. Kept a plain list otherwise
    ## (single-group output is unchanged).
    class(out) <- c("RAM_multigroup", class(out))
  }

  out
}

## Build a two-level ("within"/"between") RAM object from a lavaanify()
## parameter table that used "level: 1"/"level: 2" syntax (recognisable by
## its "level" column). See dev/PLAN-two-level-multigroup-sem.md section 4
## for the design rationale and dev/debug-two-level-multigroup/03-04 for the
## OpenMx relational (primaryKey/joinKey) translation this feeds sem().
##
## Level 1 ("within") is built exactly like the existing single-group RAM:
## its manifest variables are the real, individually-observed columns in
## the raw data, so they get an F matrix and the usual "observed variables
## free by default" mean convention -- except that convention is inverted
## here (within-level manifest means default to FIXED at 0), because the
## overall mean of a manifest variable is decomposed as
## total = between-component + within-deviation, and by that decomposition
## the within part has mean zero. Freeing it by default would duplicate,
## and be redundant with, the between-level's own free mean.
##
## Level 2 ("between") has no manifest variables of its own in this scope
## (genuine level-2-only observed covariates are not yet supported -- see
## plan). Every variable at level 2 is therefore "latent" for F-matrix
## purposes: either a "shadow" variable backing a level-1 manifest variable
## (its between-cluster component -- these default to a FREE mean, using
## the same "<var>mean" convention as the single-level default, since this
## is where a manifest variable's overall intercept lives by default), or a
## genuine level-2-only latent factor (defaults to a fixed-at-zero mean,
## the usual latent-variable convention).
##
## In both cases, these are only the DEFAULTS applied before lavaanify's
## own rows are processed: any explicit (or lavaanify-auto-added, e.g. once
## meanstructure is triggered anywhere in the model) "~1" row for a
## variable overwrites the corresponding default, exactly mirroring how the
## single-group code already lets explicit rows override its own baked-in
## defaults.
.lavaan2RAM.twolevel <- function(my.model, obs.variables, A.notation, S.notation,
                                 M.notation, A.start, S.start, M.start) {

  ## Remove constraint placeholders (.p1.==.p2. etc.) before extracting
  ## mxalgebras, exactly as the single-group code does.
  if (length(grep("^\\.", my.model$lhs)) > 0) {
    my.model <- my.model[-grep("^\\.", my.model$lhs), ]
  }

  ## Global constraints/algebras (fn1 := ..., a == b, ...) carry level==0,
  ## the same convention the single-group code uses for group==0.
  mod <- my.model[my.model$level != 0, ]

  ## Set the starting values in A/M/S if NA, across both levels together
  ## (identical logic to the single-group code).
  if (any((mod$op=="=~"|mod$op=="~")&is.na(mod$ustart))) {
    mod[(mod$op=="=~"|mod$op=="~")&is.na(mod$ustart), ]$ustart <- A.start
  }
  if (any(mod$op=="~1"&is.na(mod$ustart))) {
    mod[mod$op=="~1"&is.na(mod$ustart), ]$ustart <- M.start
  }
  if (any(mod$op=="~~"&is.na(mod$ustart)&(mod$lhs==mod$rhs))) {
    mod[mod$op=="~~"&is.na(mod$ustart)&(mod$lhs==mod$rhs), ]$ustart <- S.start
  }
  if (any(mod$op=="~~"&is.na(mod$ustart)&(mod$lhs!=mod$rhs))) {
    mod[mod$op=="~~"&is.na(mod$ustart)&(mod$lhs!=mod$rhs), ]$ustart <- 0
  }

  ## Labels: auto-generated ones get a "_within"/"_between" suffix so that,
  ## e.g., the within-level residual variance of "y1" and its between-level
  ## counterpart do not collide into a single (incorrectly equated) OpenMx
  ## parameter -- the same class of bug fixed for groups in Phase 0, but
  ## across levels instead. User-supplied labels are left untouched.
  level.suffix <- c("1"="_within", "2"="_between")
  for (i in seq_len(nrow(mod))) {
    if (mod[i, ]$label=="") {
      switch(mod[i, ]$op,
             "=~" = mod[i, ]$label <- paste0(mod[i, ]$rhs, A.notation, mod[i, ]$lhs),
             "<~" = mod[i, ]$label <- paste0(mod[i, ]$lhs, A.notation, mod[i, ]$rhs),
             "~"  = mod[i, ]$label <- paste0(mod[i, ]$lhs, A.notation, mod[i, ]$rhs),
             "~~" = mod[i, ]$label <- paste0(mod[i, ]$lhs, S.notation, mod[i, ]$rhs),
             "~1" = mod[i, ]$label <- paste0(mod[i, ]$lhs, M.notation))
      mod[i, ]$label <- paste0(mod[i, ]$label, level.suffix[[as.character(mod[i, ]$level)]])
    }
  }

  key <- with(mod, ifelse(free==0, yes=ustart, no=paste(ustart, label, sep="*")))

  ## ---------------------------------------------------------------
  ## Level 1 ("within")
  ## ---------------------------------------------------------------
  mod1 <- mod[mod$level==1, ]
  key1 <- key[mod$level==1]

  all.var1 <- unique(c(mod1$lhs, mod1$rhs))
  all.var1 <- all.var1[all.var1 != ""]
  latent1 <- unique(mod1[mod1$op=="=~", ]$lhs)
  observed1 <- all.var1[!(all.var1 %in% latent1)]

  if (!is.null(obs.variables)) {
    if (!identical(sort(observed1), sort(obs.variables))) {
      stop("Names in \"obs.variables\" do not agree with the level-1 ",
          "(within) observed variables in the model.\n")
    }
    observed1 <- obs.variables
  }

  all.var1 <- if (length(latent1) > 0) c(observed1, latent1) else observed1
  no.lat1 <- length(latent1)
  no.obs1 <- length(observed1)
  no.all1 <- no.lat1 + no.obs1

  Amatrix1 <- matrix(0, ncol=no.all1, nrow=no.all1, dimnames=list(all.var1, all.var1))
  Smatrix1 <- matrix(0, ncol=no.all1, nrow=no.all1, dimnames=list(all.var1, all.var1))
  ## Within-level variables default to a fixed-at-zero mean -- see the
  ## function-level comment above for why.
  Mmatrix1 <- matrix("0", nrow=1, ncol=no.all1, dimnames=list(1, all.var1))

  for (i in seq_len(nrow(mod1))) {
    my.line <- mod1[i, ]
    switch(my.line$op,
           "=~" = Amatrix1[my.line$rhs, my.line$lhs] <- key1[i],
           "<~" = Amatrix1[my.line$lhs, my.line$rhs] <- key1[i],
           "~"  = Amatrix1[my.line$lhs, my.line$rhs] <- key1[i],
           "~~" = Smatrix1[my.line$lhs, my.line$rhs] <-
             Smatrix1[my.line$rhs, my.line$lhs] <- key1[i],
           "~1" = Mmatrix1[1, my.line$lhs] <- key1[i]
           )
  }

  Fmatrix1 <- create.Fmatrix(c(rep(1, no.obs1), rep(0, no.lat1)), as.mxMatrix=FALSE)
  dimnames(Fmatrix1) <- list(observed1, all.var1)

  ## ---------------------------------------------------------------
  ## Level 2 ("between")
  ## ---------------------------------------------------------------
  mod2 <- mod[mod$level==2, ]
  key2 <- key[mod$level==2]

  all.var2 <- unique(c(mod2$lhs, mod2$rhs))
  all.var2 <- all.var2[all.var2 != ""]

  ## Reject genuine level-2-only observed covariates (documented as
  ## unsupported, but not actually enforced until this check was added on
  ## review): every level-2 variable must be either a "shadow" of a
  ## level-1 manifest variable, or itself a genuine latent factor (the lhs
  ## of a "=~" row at level 2, covering both a direct level-2-only
  ## covariate like "y ~ z" and a level-2-only factor indicator like
  ## "fb =~ q" where "q" never appears at level 1 either -- both would
  ## otherwise be silently absorbed as an ordinary latent variable, with
  ## no connection to any real data column, at all: no error, no warning,
  ## just a free parameter estimated in a vacuum.
  latent2 <- unique(mod2[mod2$op=="=~", ]$lhs)
  unsupported2 <- setdiff(all.var2, c(observed1, latent2))
  if (length(unsupported2) > 0) {
    stop("Genuine level-2-only observed covariates are not yet supported ",
        "by lavaan2RAM(): ", paste(unsupported2, collapse=", "), ". ",
        "Every level-2 variable must either also be a level-1 manifest ",
        "variable (its between-cluster component) or a level-2 latent ",
        "factor declared via \"=~\".\n")
  }

  no.all2 <- length(all.var2)

  Amatrix2 <- matrix(0, ncol=no.all2, nrow=no.all2, dimnames=list(all.var2, all.var2))
  Smatrix2 <- matrix(0, ncol=no.all2, nrow=no.all2, dimnames=list(all.var2, all.var2))

  ## "Shadow" variables (between-level components backing a level-1
  ## manifest variable) default to a free mean; everything else (genuine
  ## level-2-only latent factors) defaults to fixed-at-zero -- see the
  ## function-level comment above.
  is.shadow2 <- all.var2 %in% observed1
  Mmatrix2 <- matrix("0", nrow=1, ncol=no.all2, dimnames=list(1, all.var2))
  if (any(is.shadow2)) {
    Mmatrix2[1, is.shadow2] <- paste0(M.start, "*", all.var2[is.shadow2],
                                      M.notation, "_between")
  }

  for (i in seq_len(nrow(mod2))) {
    my.line <- mod2[i, ]
    switch(my.line$op,
           "=~" = Amatrix2[my.line$rhs, my.line$lhs] <- key2[i],
           "<~" = Amatrix2[my.line$lhs, my.line$rhs] <- key2[i],
           "~"  = Amatrix2[my.line$lhs, my.line$rhs] <- key2[i],
           "~~" = Smatrix2[my.line$lhs, my.line$rhs] <-
             Smatrix2[my.line$rhs, my.line$lhs] <- key2[i],
           "~1" = Mmatrix2[1, my.line$lhs] <- key2[i]
           )
  }

  ## No manifest variables at the between level in this scope (see above).
  ## Built directly (not via create.Fmatrix(rep(0, no.all2), ...)): base R's
  ## diag(0) returns an empty 0x0 matrix rather than a 1x1 zero matrix, which
  ## breaks create.Fmatrix()'s internal Diag() call whenever no.all2==1.
  Fmatrix2 <- matrix(0, nrow=0, ncol=no.all2, dimnames=list(NULL, all.var2))

  out <- list(within=list(A=Amatrix1, S=Smatrix1, F=Fmatrix1, M=Mmatrix1),
             between=list(A=Amatrix2, S=Smatrix2, F=Fmatrix2, M=Mmatrix2))

  ## Global constraints/algebras (fn1 := ..., a == b, ...), attached at the
  ## top level -- unlike the multi-group case, a two-level RAM already has
  ## a natural top-level container (it is not itself a list of "groups"),
  ## so there is no need to nest this under "within" or "between".
  if (any(my.model$level==0)) {
    mxalgebra <- list()
    con_index <- 1

    y <- my.model[my.model$level==0, , drop=FALSE]

    for (i in seq_len(nrow(y))) {
      switch(y[i, 'op'],
             ":=" = { eval(parse(text=paste0(y[i,'lhs'],
                                             "<- mxAlgebra(", y[i,'rhs'],
                                             ", name=\"", y[i,'lhs'], "\")")))
                                             eval(parse(text=paste0("mxalgebra <-
                                                  c(mxalgebra, ", y[i, 'lhs'],
                                                  "=", y[i, 'lhs'], ")")))  },
             if (y[i, 'op'] %in% c("==", ">", "<")) {
               eval(parse(text=paste0("constraint", con_index, " <- mxConstraint(", y[i, 'lhs'],
                                      y[i, 'op'], y[i, 'rhs'], ", name=\"constraint", con_index, "\")")))
               eval(parse(text=paste0("mxalgebra <- c(mxalgebra, constraint",
                                      con_index, "=constraint", con_index, ")")))
               con_index <- con_index + 1
             }
             )
    }

    out$mxalgebras <- mxalgebra
  }

  class(out) <- c("RAM_twolevel", class(out))
  out
}
