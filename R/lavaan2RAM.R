#' Convert \code{lavaan} models to RAM models
#' 
#' It converts models specified in \code{lavaan} model syntax to RAM models.
#' 
#' It uses the \code{\link[lavaan]{model.syntax}} to do the conversion.
#' Variable names and ordinary parameter labels must not contain periods.
#' Rename consistently throughout the model and data, e.g., \code{x.1} to
#' \code{x_1}, and check that the new names are unique. Definition-variable
#' references such as \code{data.x_1} remain supported: only the column name
#' following \code{data.} must be period-free. Group and cluster values may
#' contain periods.
#' Experimental: functions of parameters (`:=` in lavaan) and constraints
#' (`==`, `>`, and `<` in lavaan) will be converted to mxAlgebra and
#' mxConstraint in OpenMx. As there are differences between lavaan and OpenMx,
#' they may not work properly.
#'
#' An inequality between a parameter and a constant, e.g., \code{b > 0},
#' \code{2 < c}, or \code{lower(0)*x}, is converted to a bound on that
#' parameter rather than to an mxConstraint. This includes parameters introduced
#' only in equality expressions. The bounds are stored in the
#' \code{"lbound"} and \code{"ubound"} attributes of the output (named lists
#' keyed by parameter label) and are applied by \code{\link{sem}},
#' \code{\link{wls}} (and hence \code{\link{tssem2}}), \code{\link{osmasem}},
#' and \code{\link{osmasem2}}. \code{\link{osmasem}} ignores, with a
#' warning, bounds on parameters that it computes rather than estimates
#' (e.g., error variances of dependent variables), and its \code{A.lbound}
#' and \code{A.ubound} arguments take precedence on the \code{A} matrix.
#' Other inequalities, e.g., \code{b > 2*a}, remain
#' mxConstraints, which OpenMx satisfies within its feasibility tolerance
#' (see \code{\link[OpenMx]{mxOption}}).
#' 
#' Automatically generated observed intercepts are free by default. For a
#' free latent mean with one indicator, the indicator intercept defaults to
#' zero for compatibility. With several indicators, specify an identifying
#' intercept constraint; a warning is given when all are free. For two-level
#' models the warning is checked separately at each level and does not impose
#' an identifying constraint.
#'
#' Automatically generated parameter labels that collide with another
#' parameter are made unique with an \code{_auto} suffix and a warning.
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
#' \code{FALSE}, no marker loading is fixed automatically: fix a loading or
#' variance in the model, or set \code{auto.fix.first=TRUE}. A warning is
#' issued if no fixed nonzero loading or variance, or label constrained to a
#' nonzero constant, is found for a factor; this heuristic is not a complete
#' identification check. \code{std.lv=TRUE} fixes a latent variance at 1
#' unless the model explicitly fixes it to another value, and can override
#' a free \code{"~~"} specification in \code{model} (e.g. \code{"f ~~ c(vf1,vf2)*f"} to
#' free a factor's variance per group) with no error -- the label is simply
#' discarded, since the parameter is fixed rather than estimated. Pass
#' \code{std.lv=FALSE} whenever \code{model} identifies a latent variable's
#' scale by fixing one of its loadings instead (e.g. \code{"f =~ 1*y1"} or
#' \code{"f =~ c(1,1)*y1"} for a multiple-group model), which is
#' incompatible with also fixing its variance. \code{lavaan2RAM()} warns
#' when this happens (any parameter given an explicit label in \code{model}
#' that ends up fixed rather than free), but the warning is general --
#' checking \code{std.lv} first is faster than working backward from it.
#' @param ngroups Positive integer. The number of groups in the \code{model}.
#' A numeric character value is accepted with a deprecation warning. See
#' \code{\link[lavaan]{model.syntax}}. Parameters without an explicit label
#' in \code{model} are free and distinct in each group by default (as in
#' \code{lavaan}); use \code{lavaan}'s \code{c(label, label)} syntax to
#' constrain a parameter to be equal across groups. Has no effect (and is
#' rejected with an error) on two-level (\code{level:}) models unless
#' \code{group:} and \code{level:} blocks are nested -- see below. Groups
#' are numbered in the order of their \code{group:} blocks, whatever their
#' labels (e.g., \code{group: B} before \code{group: A}), because
#' \code{lavaan} and \code{\link{sem}} pair them with the groups in the data
#' by position. The block labels are stored in the
#' \code{"group.labels"} attribute.
#' @param auto.cov.lv.x,auto.cov.x,auto.cov.y Logical, default \code{FALSE}.
#' Request lavaan's automatic covariances for exogenous latent variables,
#' exogenous observed variables, and dependent variables, respectively.
#' @param auto.fix.first Logical, default \code{FALSE}. Fix the first loading of
#' each factor to one. With \code{std.lv=TRUE} this can impose an additional
#' restriction because the latent variance is also fixed.
#' @param auto.fix.single Logical, default \code{FALSE}. Fix the residual
#' variance of a single indicator to zero.
#' @param constraints Optional lavaan constraints using explicit parameter
#' labels. Defined parameters and equality and inequality constraints may also
#' be included in \code{model}; see the limitations above.
#' @param \dots Reserved for detecting unsupported arguments. Arbitrary
#' arguments are not passed to \code{lavaanify}; use exact argument names.
#' In particular, \code{orthogonal}, \code{warn}, \code{debug},
#' \code{effect.coding}, \code{group.equal}, and \code{correlation} are
#' not supported. Specify shared group parameters with explicit labels.
#' @param auto Reserved to detect an unsupported argument; do not supply it.
#' It is not lavaan's \code{auto} option.
#' @details Continuous RAM syntax supports regressions, loadings, variances,
#' covariances, means, and the constraints described above. Ordinal thresholds,
#' composites, exploratory factor blocks, scaling operators, and group weights
#' are not supported and are rejected. The returned RAM describes a covariance
#' model; correlation restrictions are handled by the fitting function (for
#' example, \code{cor.analysis} in \code{wls} or \code{osmasem2}).
#' @return A list of RAM specification with \code{A}, \code{S}, \code{F}, and
#' \code{M} matrices. When \code{ngroups > 1}, a list of such lists (one per
#' group) is returned instead, with class \code{"RAM_multigroup"}. When
#' \code{model} uses \code{lavaan}'s two-level syntax (\code{level: 1} /
#' \code{level: 2} blocks), a list with \code{within} and \code{between}
#' elements (each themselves a list of \code{A}, \code{S}, \code{F}, and
#' \code{M} matrices) is returned instead, with class \code{"RAM_twolevel"}.
#' Level-1 (within) manifest variables that also appear at level 2 default
#' to a fixed-at-zero mean (the overall intercept lives at the between level
#' by convention), whereas within-only observed variables, which have no
#' between component, default to a free level-1 mean, as in \code{lavaan};
#' level-2
#' (between) "shadow" variables backing a level-1 manifest variable default
#' to a free mean, while genuine level-2-only latent factors default to
#' fixed-at-zero, mirroring the usual latent-variable convention. Genuine
#' level-2-only observed covariates (e.g. a real cluster-level covariate
#' with its own data column, as opposed to a level-1 manifest variable's
#' between-cluster component) are supported when their data values are
#' constant within each cluster. Models with nested \code{group:} and
#' \code{level:} blocks and \code{ngroups > 1} return a list of two-level
#' RAM objects with class \code{"RAM_multigroup_twolevel"}.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[semPlot]{ramModel}}, \code{\link[metaSEM]{Becker92}},
#' \code{\link[metaSEM]{Becker09}}, \code{\link[metaSEM]{Digman97}},
#' \code{\link[metaSEM]{Hunter83}}, \code{\link[metaSEM]{as.mxMatrix}},
#' \code{\link[metaSEM]{checkRAM}}
#' @keywords methods tssem
#' @examples
#' ## Request an automatic residual covariance between dependent variables.
#' lavaan2RAM("y ~ x\nz ~ x", auto.cov.y=TRUE)
#' ## Choose a marker-loading scale when the latent variance is not fixed.
#' lavaan2RAM("f =~ y1 + y2 + y3", std.lv=FALSE, auto.fix.first=TRUE)
#' ## Parameter constraints use explicit labels.
#' lavaan2RAM("y ~ b*x", constraints="b == 0.5")
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
                       ngroups=1, auto.cov.lv.x=FALSE, auto.cov.x=FALSE,
                       auto.cov.y=FALSE, auto.fix.first=FALSE,
                       auto.fix.single=FALSE, constraints=NULL, ..., auto=NULL) {
  ## if (!requireNamespace("lavaan", quietly=TRUE))    
  ##     stop("\"lavaan\" package is required for this function.")

  ## Inspect names before forcing arguments, so an abbreviated argument cannot
  ## silently select a different option through R's partial matching.
  raw.args <- as.list(sys.call())[-1L]
  raw.names <- names(raw.args)
  ## A literal ... in this call is expanded by R before matching arguments,
  ## but sys.call() still contains only the symbol. Inspect the caller's
  ## unevaluated dots to retain their original names, also through nested
  ## forwarding wrappers and without evaluating their values.
  forwarded <- vapply(seq_along(raw.args), function(i)
    identical(unname(raw.args[i]), list(as.name("..."))), logical(1))
  if (any(forwarded)) {
    dots <- eval(quote(base::substitute(list(...))), envir=parent.frame())
    dot.names <- names(as.list(dots)[-1L])
    raw.names <- c(raw.names, dot.names)
  }
  if (length(raw.names)) {
    raw.names <- raw.names[nzchar(raw.names)]
    allowed <- names(formals(lavaan2RAM))
    partial <- setdiff(raw.names, allowed)
    if (length(partial))
      stop("Unsupported argument(s): ", paste(unique(partial), collapse=", "),
           ". Use the exact documented argument names.", call.=FALSE)
  }
  if (!missing(auto))
    stop("'auto' is not supported; specify the documented automatic options individually.", call.=FALSE)
  extra <- as.list(substitute(list(...)))[-1L]
  if (length(extra)) {
    extra.names <- names(extra)
    if (is.null(extra.names)) extra.names <- rep("", length(extra))
    extra.names[!nzchar(extra.names)] <- "<unnamed>"
    stop("Unsupported argument(s): ", paste(unique(extra.names), collapse=", "),
         ". See ?lavaan2RAM for supported options.", call.=FALSE)
  }
  flags <- list(auto.var=auto.var, std.lv=std.lv, auto.cov.lv.x=auto.cov.lv.x,
                auto.cov.x=auto.cov.x, auto.cov.y=auto.cov.y,
                auto.fix.first=auto.fix.first, auto.fix.single=auto.fix.single)
  valid.flag <- vapply(flags, function(x) is.logical(x) && length(x)==1L && !is.na(x), logical(1))
  if (any(!valid.flag))
    stop("These arguments must be TRUE or FALSE: ",
         paste(names(flags)[!valid.flag], collapse=", "), call.=FALSE)

  if (is.character(ngroups) && length(ngroups) == 1L && !is.na(ngroups)) {
    numeric.ngroups <- suppressWarnings(as.numeric(ngroups))
    if (is.finite(numeric.ngroups) && numeric.ngroups >= 1 &&
        numeric.ngroups == floor(numeric.ngroups) && numeric.ngroups <= .Machine$integer.max) {
      warning("Character 'ngroups' is deprecated; supply a positive integer instead.", call.=FALSE)
      ngroups <- numeric.ngroups
    }
  }
  if (!is.numeric(ngroups) || is.complex(ngroups) || length(ngroups) != 1L ||
      is.na(ngroups) || !is.finite(ngroups) || ngroups < 1 ||
      ngroups != floor(ngroups) || ngroups > .Machine$integer.max)
    stop("'ngroups' must be a single positive integer.", call.=FALSE)
  ngroups <- as.integer(ngroups)

  ## Default: fix the latent independent variables at 1
  my.model <- lavaan::lavaanify(model, fixed.x = FALSE, auto.var=auto.var,
                                std.lv=std.lv, ngroups=ngroups,
                                auto.cov.lv.x=auto.cov.lv.x, auto.cov.x=auto.cov.x,
                                auto.cov.y=auto.cov.y, auto.fix.first=auto.fix.first,
                                auto.fix.single=auto.fix.single, constraints=constraints)

  unsupported <- setdiff(unique(my.model$op), c("=~", "~", "~~", "~1", ":=", "==", "<", ">"))
  if (length(unsupported))
    stop("Unsupported lavaan operator(s): ", paste(unsupported, collapse=", "),
         ". Ordinal thresholds, scaling, composites, and group weights are not supported.", call.=FALSE)
  if ("efa" %in% names(my.model) && any(nzchar(my.model$efa), na.rm=TRUE))
    stop("Exploratory factor analysis blocks (efa) are not supported.", call.=FALSE)

  paths <- my.model$op %in% c("=~", "~", "~~", "~1")
  .check.no.period(c(my.model$lhs[paths], my.model$rhs[paths], obs.variables),
                   "Model variable names")
  .check.parameter.names(my.model$label[paths | my.model$op == ":="])
  .check.no.period(c(A.notation, S.notation, M.notation), "Parameter notation")

  ## lavaan adds placeholder equalities for parameters already sharing an
  ## explicit label. Remove only equalities proven redundant; a placeholder
  ## in any other constraint must not be silently discarded.
  placeholder <- grepl("^\\.p[0-9]+\\.$", my.model$lhs) & my.model$op == "=="
  for (i in which(placeholder)) {
    left <- which(my.model$plabel == my.model$lhs[i] & my.model$op != "==")
    right <- which(my.model$plabel == my.model$rhs[i] & my.model$op != "==")
    redundant <- length(left)==1L && length(right)==1L &&
      nzchar(my.model$label[left]) && my.model$label[left] == my.model$label[right]
    if (!redundant)
      stop("Unsupported generated equality constraint: ", my.model$lhs[i],
           " == ", my.model$rhs[i], ". Use explicit named parameters.", call.=FALSE)
  }
  if (any(placeholder)) my.model <- my.model[!placeholder, , drop=FALSE]
  constraint.rows <- my.model$op %in% c(":=", "==", "<", ">")
  if (any(grepl("\\.p[0-9]+\\.", paste(my.model$lhs[constraint.rows],
                                             my.model$rhs[constraint.rows]))))
    stop("Constraints referencing lavaan's automatic .p labels are not supported; use explicit parameter labels.", call.=FALSE)
  symbols <- character()
  if (any(constraint.rows)) {
    expressions <- c(my.model$lhs[constraint.rows], my.model$rhs[constraint.rows])
    symbols <- unlist(lapply(expressions, function(x) {
      parsed <- tryCatch(parse(text=x), error=function(e) NULL)
      if (is.null(parsed)) character() else all.vars(parsed)
    }), use.names=FALSE)
    ## Constraint expressions can contain qualified OpenMx references such
    ## as aux.value; ordinary labels on model paths were checked above.
    .check.parameter.names(symbols, free=rep(FALSE, length(symbols)))
  }

  ## lavaan pairs group blocks with the groups in the data by position,
  ## whatever their labels (e.g., "group: B" before "group: A", or "group: 2"
  ## before "group: 1"). Number the groups by their order in the syntax; rows
  ## shared by all groups (e.g., constraints) get 0. lavaanify() gives text
  ## labels for "group: A" syntax, with "" for the shared rows. The labels of
  ## "group:" blocks are kept for sem() to compare with the data.
  group.labels <- NULL
  if (!is.null(my.model$group)) {
    group <- as.character(my.model$group)
    shared <- is.na(group) | group %in% c("", "0")
    labels <- unique(group[!shared])
    my.model$group <- ifelse(shared, 0L, match(group, labels))
    if (any(grepl("^[[:space:]]*group[[:space:]]*:",
                  unlist(strsplit(model, "[\n;]")))))
      group.labels <- labels
  }

  if (!std.lv) {
    ## Explicit label == nonzero-constant constraints also anchor a scale.
    ## Do not infer identification from general expressions or equalities
    ## between two unknown parameters.
    equalities <- my.model[my.model$op == "==", , drop=FALSE]
    left.number <- suppressWarnings(as.numeric(equalities$lhs))
    right.number <- suppressWarnings(as.numeric(equalities$rhs))
    anchored.labels <- c(equalities$lhs[is.finite(right.number) & right.number != 0],
                         equalities$rhs[is.finite(left.number) & left.number != 0])
    anchored.labels <- anchored.labels[nzchar(anchored.labels)]
    loading.rows <- my.model[my.model$op == "=~", , drop=FALSE]
    if (nrow(loading.rows)) {
      group.key <- if ("group" %in% names(my.model)) as.character(my.model$group) else rep("1", nrow(my.model))
      level.key <- if ("level" %in% names(my.model)) as.character(my.model$level) else rep("1", nrow(my.model))
      loading.index <- which(my.model$op == "=~")
      factors <- unique(data.frame(name=my.model$lhs[loading.index],
                                   group=group.key[loading.index], level=level.key[loading.index],
                                   stringsAsFactors=FALSE))
      unanchored <- character()
      for (i in seq_len(nrow(factors))) {
        factor <- factors$name[i]
        same <- group.key == factors$group[i] & level.key == factors$level[i]
        loading <- same & my.model$op == "=~" & my.model$lhs == factor
        variance <- same & my.model$op == "~~" & my.model$lhs == factor & my.model$rhs == factor
        fixed.nonzero <- my.model$free == 0 & !is.na(my.model$ustart) & my.model$ustart != 0
        constrained.nonzero <- nzchar(my.model$label) & my.model$label %in% anchored.labels
        if (!any((fixed.nonzero | constrained.nonzero) & (loading | variance)))
          unanchored <- c(unanchored, factor)
      }
      if (length(unanchored))
        warning("With std.lv=FALSE, factor(s) ", paste(unique(unanchored), collapse=", "),
                " have no fixed nonzero loading or variance in at least one group/level; ",
                "their scales may be unidentified. Specify an identifying constraint or auto.fix.first=TRUE.",
                call.=FALSE)
    }
  }

  ## Warn about a parameter that the user gave an explicit label but that
  ## lavaanify() has fixed anyway (free==0): a label on a fixed value does
  ## nothing (there is no free parameter for it to equate or report), so
  ## this is virtually always unintentional. By far the most common cause
  ## is a latent variable's own "~~" variance/covariance clashing with
  ## std.lv=TRUE (this function's default), which fixes a labelled free
  ## latent variance at 1 for identification, silently discarding that
  ## explicit "~~" label with no error -- e.g. writing
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
           "(the default here), which can fix a labelled free variance at 1 -- ",
           "pass std.lv=FALSE if you ",
           "intended it to be free (e.g. because the model already fixes a ",
           "loading for identification instead).", call.=FALSE)
  }

  ## Inequalities between a free parameter and a constant (e.g. "b > 3") are
  ## bounds. lavaanify() moves some of them into its "lower"/"upper" columns
  ## and keeps others as "<"/">" rows (e.g. in multiple-group models); both
  ## are recorded on the RAM object so that sem() can apply them.
  bounds <- .lavaan.bounds(my.model)
  my.model <- bounds$model

  ## Reserve explicit labels and defined-algebra names before generating any
  ## labels. The same key lets a default mean and its explicit lavaan row
  ## refer to the same generated parameter.
  explicit <- unique(c(my.model$label[nzchar(my.model$label)],
                       my.model$lhs[my.model$op == ":="]))
  used <- explicit
  assigned <- new.env(parent=emptyenv())
  labeler <- function(candidate, key) {
    if (exists(key, envir=assigned, inherits=FALSE)) return(get(key, envir=assigned))
    result <- candidate
    if (result %in% used) {
      if (candidate %in% symbols && !(candidate %in% explicit))
        stop("Generated label '", candidate, "' is ambiguous in a constraint; assign explicit distinct labels.", call.=FALSE)
      k <- 1L
      repeat {
        result <- paste0(candidate, "_auto", k)
        if (!(result %in% used)) break
        k <- k+1L
      }
      warning("Generated parameter label '", candidate, "' conflicts with another label; using '", result, "'.", call.=FALSE)
    }
    used <<- c(used, result)
    assign(key, result, envir=assigned)
    result
  }

  ## Two-level ("level: 1" / "level: 2") syntax: lavaanify() reports this
  ## via a "level" column. A "group" column is present as well only with
  ## nested "group:"/"level:" blocks, which give one two-level RAM per group.
  if ("level" %in% colnames(my.model)) {
    if ("group" %in% colnames(my.model)) {
      groups <- sort(unique(my.model$group[my.model$group > 0]))
      if (length(groups) != ngroups)
        stop("The number of group blocks does not match 'ngroups'.")
      out <- lapply(groups, function(g) {
        rows <- my.model$group == g | my.model$group == 0
        .lavaan2RAM.twolevel(my.model[rows, , drop=FALSE],
                             obs.variables=obs.variables,
                             A.notation=A.notation, S.notation=S.notation,
                             M.notation=M.notation, A.start=A.start,
                             S.start=S.start, M.start=M.start,
                             group.suffix=paste0("_g", g), labeler=labeler)
      })
      global.algebras <- out[[1]]$mxalgebras
      for (i in seq_along(out)) out[[i]]$mxalgebras <- NULL
      attr(out, "mxalgebras") <- global.algebras
      class(out) <- c("RAM_multigroup_twolevel", class(out))
      return(.attach.bounds(out, bounds, group.labels))
    }
    if (ngroups > 1) {
      stop("Combined group and level models require nested 'group:' and ",
           "'level:' blocks in the lavaan syntax.")
    }
    return(.attach.bounds(
      .lavaan2RAM.twolevel(my.model, obs.variables=obs.variables,
                           A.notation=A.notation, S.notation=S.notation,
                           M.notation=M.notation, A.start=A.start,
                           S.start=S.start, M.start=M.start, labeler=labeler), bounds))
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
    ## Prepare the labels
    for (i in seq_len(nrow(mod))) {
      ## if there is no label
      if (mod[i, ]$label=="") {
        switch(mod[i, ]$op,
               "=~" = mod[i, ]$label <- paste0(mod[i, ]$rhs, A.notation, mod[i, ]$lhs),
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
        if (mod[i, ]$free != 0) {
          key.id <- paste(gp, "single", mod[i, ]$op, mod[i, ]$lhs, mod[i, ]$rhs, sep="|")
          mod[i, ]$label <- labeler(mod[i, ]$label, key.id)
        }
      }
    }

    mean.labels <- vapply(observed, function(v) {
      overridden <- any(mod$op == "~1" & mod$lhs == v & (mod$user != 0 | mod$free != 0))
      if (overridden) return("0")
      paste0(M.start, "*", labeler(paste0(v, M.notation, mean.suffix),
          paste(gp, "single", "~1", v, "", sep="|")))
    }, character(1))
    Mmatrix <- matrix(c(mean.labels, rep("0", no.lat)),
                      nrow=1, ncol=no.all, dimnames=list(1, all.var))

    ## Replace NA to 0 in ustart if there are still NA
    ## mod$ustart[is.na(mod$ustart)] <- 0
        
    ## keys in as.mxMatrix format
    key <- ifelse(mod$free == 0, yes=mod$ustart,
                  no=paste(mod$ustart, mod$label, sep="*"))
        
    ## Preserve the established single-indicator marker convention when a
    ## latent mean is free. For other indicators, lavaanify's automatic
    ## fixed-zero intercept rows must not override the free defaults.
    free.latent.means <- mod$lhs[mod$op == "~1" & mod$lhs %in% latent & mod$free > 0]
    marker.indicators <- unique(unlist(lapply(free.latent.means, function(factor) {
      indicators <- unique(mod$rhs[mod$op == "=~" & mod$lhs == factor])
      if (length(indicators) == 1L) indicators else character()
    })))
    for (i in seq_len(nrow(mod))) {
      my.line <- mod[i, ]
      if (my.line$op == "~1" && my.line$user == 0 && my.line$free == 0 &&
          my.line$lhs %in% observed && !my.line$lhs %in% marker.indicators) next
      switch(my.line$op,
             ## lhs: IV; rhs: DV
             "=~" = Amatrix[my.line$rhs, my.line$lhs] <- key[i],
             ## lhs: DV; rhs: IV
             "~"  = Amatrix[my.line$lhs, my.line$rhs] <- key[i],
             "~~" = Smatrix[my.line$lhs, my.line$rhs] <-
               Smatrix[my.line$rhs, my.line$lhs] <- key[i],
             ## means
             "~1" = Mmatrix[1, my.line$lhs] <- key[i]
             )  ## from switch
    }  ## from for loop

    ## An unrestricted latent mean and unrestricted indicator intercepts
    ## may be redundant. Do not silently impose further restrictions.
    mean.free <- setNames(as.mxMatrix(Mmatrix, name="checkMeans")$free[1, ],
                          colnames(Mmatrix))
    for (factor in unique(free.latent.means)) {
      indicators <- unique(mod$rhs[mod$op == "=~" & mod$lhs == factor])
      if (length(indicators) > 0L && all(mean.free[indicators]))
        warning("The mean of factor '", factor, "' and all its indicator intercepts are free; ",
                "the mean structure may be unidentified. Fix an indicator intercept or provide an identifying constraint.", call.=FALSE)
    }

    Fmatrix <- create.Fmatrix(c(rep(1, no.obs), rep(0, no.lat)), as.mxMatrix=FALSE)
    dimnames(Fmatrix) <- list(observed, all.var)
    
    out[[gp]] <- list(A=Amatrix, S=Smatrix, F=Fmatrix, M=Mmatrix)
  }

  ## Add group names, 1, 2, 3... to the list
  names(out) <- seq_along(out)
    
  ## Check if there are constraints or algebras
  if (any(my.model$group==0)) {
        
    ## Constraints and algebras only
    mxalgebra <- .lavaan.algebras(my.model[my.model$group==0, , drop=FALSE],
                                  my.model)

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

  .attach.bounds(out, bounds, group.labels)
}

## Bounds implied by lavaan syntax. Inequality rows comparing a free
## parameter's label with a numeric constant are converted to bounds and
## removed from the parameter table; any other "<"/">" row (e.g. "b > 2*a")
## remains an mxConstraint. Free parameters that lavaanify() has already
## bounded through its "lower"/"upper" columns, which may be unlabeled (e.g.
## "lower(0)*x"), are returned as rows to be resolved to their RAM labels by
## .attach.bounds().
.lavaan.bounds <- function(my.model) {
  lbound <- ubound <- numeric(0)
  add <- function(bounds, label, value, keep) {
    if (label %in% names(bounds)) value <- keep(bounds[[label]], value)
    bounds[label] <- value
    bounds
  }
  ## Evaluate a constant such as "0" or "-1 + 0.5" with arithmetic only.
  arithmetic <- list2env(list("+"=`+`, "-"=`-`, "*"=`*`, "/"=`/`, "^"=`^`,
                              "("=`(`, exp=exp, log=log, sqrt=sqrt),
                         parent=emptyenv())
  constant <- function(x) {
    value <- tryCatch(eval(parse(text=x)[[1]], arithmetic),
                      error=function(e) NULL)
    if (is.numeric(value) && length(value) == 1L && is.finite(value)) value
    else NULL
  }
  labels <- unique(my.model$label[my.model$free > 0 & nzchar(my.model$label)])
  ## Equality substitutions can introduce parameters absent from A/S/M.
  equality.rows <- my.model$op == "=="
  introduced <- all.vars(parse(text=c(my.model$lhs[equality.rows],
                                      my.model$rhs[equality.rows])))
  labels <- union(labels, introduced[!startsWith(introduced, "data.")])
  converted <- logical(nrow(my.model))
  for (i in which(my.model$op %in% c("<", ">"))) {
    lhs <- my.model$lhs[i]
    rhs <- my.model$rhs[i]
    greater <- my.model$op[i] == ">"
    if (lhs %in% labels && !is.null(value <- constant(rhs))) {
      label <- lhs
    } else if (rhs %in% labels && !is.null(value <- constant(lhs))) {
      label <- rhs
      greater <- !greater
    } else next
    if (greater) lbound <- add(lbound, label, value, max)
    else ubound <- add(ubound, label, value, min)
    converted[i] <- TRUE
  }
  ## lavaanify() adds a "lower" and/or an "upper" column only when needed.
  for (column in c("lower", "upper"))
    if (is.null(my.model[[column]])) my.model[[column]] <- NA_real_
  rows <- my.model[my.model$free > 0 & (is.finite(my.model$lower) |
                                          is.finite(my.model$upper)), ,
                   drop=FALSE]
  ## Position of each row's group in the RAM list. The "group" column holds
  ## group labels for "group: A" syntax; lavaan2RAM() orders groups the same.
  rows$group.index <- if (is.null(my.model$group)) rep(NA_integer_, nrow(rows))
    else match(rows$group, sort(unique(my.model$group[my.model$group > 0])))
  list(model=my.model[!converted, , drop=FALSE], lbound=lbound,
       ubound=ubound, rows=rows)
}

## Attach the bounds from .lavaan.bounds() to a RAM object as "lbound" and
## "ubound" attributes (named lists keyed by parameter label, as in sem()).
## Bounded rows are resolved to the label written in their RAM matrix cell.
## The labels of "group:" blocks of a multiple-group RAM are attached as
## "group.labels".
.attach.bounds <- function(out, bounds, group.labels=NULL) {
  if (!is.null(group.labels) &&
      inherits(out, c("RAM_multigroup", "RAM_multigroup_twolevel")))
    attr(out, "group.labels") <- group.labels
  lbound <- bounds$lbound
  ubound <- bounds$ubound
  rows <- bounds$rows
  for (i in seq_len(nrow(rows))) {
    ram <- out
    if (inherits(out, c("RAM_multigroup", "RAM_multigroup_twolevel")))
      ram <- out[[rows$group.index[i]]]
    if (!is.null(rows$level) && is.list(ram$within))
      ram <- ram[[if (rows$level[i] == 1) "within" else "between"]]
    lhs <- rows$lhs[i]
    rhs <- rows$rhs[i]
    cell <- tryCatch(switch(rows$op[i],
                            "~"=ram$A[lhs, rhs], "=~"=ram$A[rhs, lhs],
                            "~~"=ram$S[lhs, rhs], "~1"=ram$M[1, lhs]),
                     error=function(e) NULL)
    if (is.null(cell) || !grepl("*", cell, fixed=TRUE)) {
      warning("The bound on '", lhs, " ", rows$op[i], " ", rhs,
              "' could not be matched to a RAM parameter and is ignored.",
              call.=FALSE)
      next
    }
    label <- sub("^[^*]*\\*", "", cell)
    if (is.finite(rows$lower[i]))
      lbound[label] <- max(c(lbound[label], rows$lower[i]), na.rm=TRUE)
    if (is.finite(rows$upper[i]))
      ubound[label] <- min(c(ubound[label], rows$upper[i]), na.rm=TRUE)
  }
  if (length(lbound)) attr(out, "lbound") <- as.list(lbound)
  if (length(ubound)) attr(out, "ubound") <- as.list(ubound)
  out
}

## mxAlgebra (":=") and mxConstraint ("==", "<", ">") objects from the rows
## of a lavaanify() parameter table that are shared by all groups or levels,
## as a list named by the defined parameters and "constraint1",
## "constraint2", ... "model" is the full parameter table, used to reject
## defined-parameter names that would clash with other model entities.
.lavaan.algebras <- function(rows, model) {
  defined <- rows$lhs[rows$op == ":="]
  reserved <- c("Amatrix", "Smatrix", "Fmatrix", "Mmatrix", "Id", "Id_A",
                "expCov", "expMean", "data", "fitfunction", "expectation",
                "compute")
  labels <- unique(model$label[nzchar(model$label) & model$op != ":="])
  clash <- defined[defined %in% c(reserved, labels) |
                     grepl("^constraint[0-9]+$", defined)]
  if (length(clash))
    stop("Defined parameter name(s) ", paste(unique(clash), collapse=", "),
         " clash with parameter labels, the constraint names \"constraint1\", ",
         "\"constraint2\", ..., or the names of matrices and algebras in the ",
         "fitted model. Please rename them.", call.=FALSE)
  out <- list()
  count <- 0L
  for (i in seq_len(nrow(rows))) {
    op <- rows$op[i]
    if (op == ":=") {
      out[[rows$lhs[i]]] <- mxAlgebraFromString(rows$rhs[i], name=rows$lhs[i])
    } else if (op %in% c("==", "<", ">")) {
      count <- count + 1L
      name <- paste0("constraint", count)
      ## Spaces around the operator: "c<-5*a" would parse as an assignment
      formula <- str2lang(paste(rows$lhs[i], op, rows$rhs[i]))
      out[[name]] <- eval(bquote(mxConstraint(.(formula), name=.(name))))
    }
  }
  out
}

## Build a two-level ("within"/"between") RAM object from a lavaanify()
## parameter table that used "level: 1"/"level: 2" syntax (recognisable by
## its "level" column). See dev/technical.md for the current design and
## dev/debug-two-level-multigroup/03-04 for the original
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
## and be redundant with, the between-level's own free mean. A within-only
## observed variable (one that does not appear at level 2) has no between
## component, so its mean is free at level 1, as in lavaan.
##
## Level 2 ("between") variables are either a "shadow" variable backing a
## level-1 manifest variable (its between-cluster component -- these
## default to a FREE mean, using the same "<var>mean" convention as the
## single-level default, since this is where a manifest variable's overall
## intercept lives by default), an observed cluster-level covariate (a
## level-2-only observed variable, also with a free mean, whose values must
## be constant within each cluster), or a genuine level-2-only latent
## factor (defaults to a fixed-at-zero mean, the usual latent-variable
## convention). Only the cluster-level covariates are manifest in the
## between-level F matrix.
##
## In both cases, these are only the DEFAULTS applied before lavaanify's
## own rows are processed: any explicit (or lavaanify-auto-added, e.g. once
## meanstructure is triggered anywhere in the model) "~1" row for a
## variable overwrites the corresponding default, exactly mirroring how the
## single-group code already lets explicit rows override its own baked-in
## defaults.
.lavaan2RAM.twolevel <- function(my.model, obs.variables, A.notation, S.notation,
                                 M.notation, A.start, S.start, M.start,
                                 group.suffix="", labeler) {

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
  level.suffix <- c("1"=paste0("_within", group.suffix),
                    "2"=paste0("_between", group.suffix))
  for (i in seq_len(nrow(mod))) {
    if (mod[i, ]$label=="") {
      switch(mod[i, ]$op,
             "=~" = mod[i, ]$label <- paste0(mod[i, ]$rhs, A.notation, mod[i, ]$lhs),
             "~"  = mod[i, ]$label <- paste0(mod[i, ]$lhs, A.notation, mod[i, ]$rhs),
             "~~" = mod[i, ]$label <- paste0(mod[i, ]$lhs, S.notation, mod[i, ]$rhs),
             "~1" = mod[i, ]$label <- paste0(mod[i, ]$lhs, M.notation))
      mod[i, ]$label <- paste0(mod[i, ]$label, level.suffix[[as.character(mod[i, ]$level)]])
      if (mod[i, ]$free != 0) {
        key.id <- paste(group.suffix, mod[i, ]$level, mod[i, ]$op, mod[i, ]$lhs, mod[i, ]$rhs, sep="|")
        mod[i, ]$label <- labeler(mod[i, ]$label, key.id)
      }
    }
  }

  key <- ifelse(mod$free == 0, yes=mod$ustart,
                  no=paste(mod$ustart, mod$label, sep="*"))

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
  ## Within-level variables that also have a between component default to
  ## a fixed-at-zero mean -- see the function-level comment above for why.
  ## A within-only observed variable has no between component, so its mean
  ## is free at level 1, as in lavaan.
  Mmatrix1 <- matrix("0", nrow=1, ncol=no.all1, dimnames=list(1, all.var1))
  level2.vars <- unique(c(mod$lhs[mod$level==2], mod$rhs[mod$level==2]))
  within.only <- setdiff(observed1, level2.vars)
  if (length(within.only)) {
    labels <- vapply(within.only, function(v) {
      overridden <- any(mod1$op == "~1" & mod1$lhs == v & (mod1$user != 0 | mod1$free != 0))
      if (overridden) return("0")
      paste0(M.start, "*", labeler(paste0(v, M.notation, "_within", group.suffix),
          paste(group.suffix, 1, "~1", v, "", sep="|")))
    }, character(1))
    Mmatrix1[1, within.only] <- labels
  }

  ## lavaanify() adds a fixed-at-zero "~1" row (user==0, free==0) for every
  ## variable once any "~1" appears in the model. These rows would override
  ## the default means set above (e.g., fix the between-level means of all
  ## other variables at 0), so they are skipped, as in the single-group code.
  ## Explicit "~1" rows in the syntax (user==1) still override the defaults.
  auto.mean <- function(x) x$op=="~1" && x$user==0 && x$free==0

  for (i in seq_len(nrow(mod1))) {
    my.line <- mod1[i, ]
    if (auto.mean(my.line)) next
    switch(my.line$op,
           "=~" = Amatrix1[my.line$rhs, my.line$lhs] <- key1[i],
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

  ## Variables found only at level 2 and not declared as factors are
  ## observed cluster-level columns. The fitter checks that each is constant
  ## within a cluster before putting one value per cluster in between.data.
  latent2 <- unique(mod2[mod2$op=="=~", ]$lhs)
  observed2 <- setdiff(all.var2, c(observed1, latent2))

  no.all2 <- length(all.var2)

  Amatrix2 <- matrix(0, ncol=no.all2, nrow=no.all2, dimnames=list(all.var2, all.var2))
  Smatrix2 <- matrix(0, ncol=no.all2, nrow=no.all2, dimnames=list(all.var2, all.var2))

  ## "Shadow" variables (between-level components backing a level-1
  ## manifest variable) default to a free mean; everything else (genuine
  ## level-2-only latent factors) defaults to fixed-at-zero -- see the
  ## function-level comment above.
  is.shadow2 <- all.var2 %in% c(observed1, observed2)
  Mmatrix2 <- matrix("0", nrow=1, ncol=no.all2, dimnames=list(1, all.var2))
  if (any(is.shadow2)) {
    labels <- vapply(all.var2[is.shadow2], function(v) {
      overridden <- any(mod2$op == "~1" & mod2$lhs == v & (mod2$user != 0 | mod2$free != 0))
      if (overridden) return("0")
      paste0(M.start, "*", labeler(paste0(v, M.notation, "_between", group.suffix),
          paste(group.suffix, 2, "~1", v, "", sep="|")))
    }, character(1))
    Mmatrix2[1, is.shadow2] <- labels
  }

  for (i in seq_len(nrow(mod2))) {
    my.line <- mod2[i, ]
    if (auto.mean(my.line)) next
    switch(my.line$op,
           "=~" = Amatrix2[my.line$rhs, my.line$lhs] <- key2[i],
           "~"  = Amatrix2[my.line$lhs, my.line$rhs] <- key2[i],
           "~~" = Smatrix2[my.line$lhs, my.line$rhs] <-
             Smatrix2[my.line$rhs, my.line$lhs] <- key2[i],
           "~1" = Mmatrix2[1, my.line$lhs] <- key2[i]
           )
  }

  ## A free latent mean and free intercepts for all of its indicators can
  ## represent the same location shift. Report this without changing the
  ## user's model; simple equalities fixing a mean are recognised.
  eq <- my.model[my.model$op == "==", , drop=FALSE]
  lhs.num <- suppressWarnings(as.numeric(eq$lhs))
  rhs.num <- suppressWarnings(as.numeric(eq$rhs))
  fixed.labels <- c(eq$lhs[is.finite(rhs.num)], eq$rhs[is.finite(lhs.num)])
  zero.labels <- c(eq$lhs[is.finite(rhs.num) & rhs.num == 0],
                   eq$rhs[is.finite(lhs.num) & lhs.num == 0])
  check.means <- function(part, Mmatrix, level) {
    mm <- as.mxMatrix(Mmatrix, name="checkLevelMeans")
    free <- setNames(mm$free[1, ], colnames(Mmatrix))
    labs <- setNames(mm$labels[1, ], colnames(Mmatrix))
    fixed <- !is.na(labs) & labs %in% fixed.labels
    free[fixed] <- FALSE
    factors <- unique(part$lhs[part$op == "=~"])
    for (factor in factors) {
      loadings <- part[part$op == "=~" & part$lhs == factor, , drop=FALSE]
      ## A known fixed-zero loading cannot anchor the factor mean. Free
      ## loadings remain possible anchors even when their starting value is 0.
      definitely.zero <- (loadings$free == 0 & !is.na(loadings$ustart) & loadings$ustart == 0) |
        (!is.na(loadings$label) & loadings$label %in% zero.labels)
      indicators <- unique(loadings$rhs[!definitely.zero])
      if (isTRUE(free[[factor]]) && (!length(indicators) || all(free[indicators])))
        warning("The ", level, "-level mean of factor '", factor,
                "' has no fixed indicator intercept with a nonzero loading; the mean structure ",
                "may be unidentified. Fix an indicator intercept or provide an identifying constraint.",
                call.=FALSE)
    }
  }
  check.means(mod1, Mmatrix1, "within")
  check.means(mod2, Mmatrix2, "between")

  Fmatrix2 <- matrix(0, nrow=length(observed2), ncol=no.all2,
                     dimnames=list(observed2, all.var2))
  if (length(observed2)) Fmatrix2[cbind(seq_along(observed2), match(observed2, all.var2))] <- 1

  out <- list(within=list(A=Amatrix1, S=Smatrix1, F=Fmatrix1, M=Mmatrix1),
             between=list(A=Amatrix2, S=Smatrix2, F=Fmatrix2, M=Mmatrix2))

  ## Global constraints/algebras (fn1 := ..., a == b, ...), attached at the
  ## top level -- unlike the multi-group case, a two-level RAM already has
  ## a natural top-level container (it is not itself a list of "groups"),
  ## so there is no need to nest this under "within" or "between".
  if (any(my.model$level==0)) {
    mxalgebra <- .lavaan.algebras(my.model[my.model$level==0, , drop=FALSE],
                                  my.model)
    out$mxalgebras <- mxalgebra
  }

  class(out) <- c("RAM_twolevel", class(out))
  out
}
