## Main function to convert a lavaan parameter table to RAM
.lavaan2RAM <- function(x, obs.variables = NULL, A.notation="ON", S.notation="WITH",
                        M.notation="mean") {
    ## replace "." with "" as OpenMx does not like "." in model
    x$plabel <- gsub("\\.", "", x$plabel)
    ## remove constraints 
    ## there are constraints such as .p1. and .p2; remove them
    if (length(grep("^\\.", x$lhs)) >0 ) x <- x[-grep("^\\.", x$lhs), ]
    ## label as the variable labels
    #my.model$label <- with(my.model, ifelse(label=="", yes=plabel, no=label))
    
    ## set the starting values in A as 0 if NA
    if (any((x$op=="=~"|x$op=="~")&is.na(x$ustart))) {
        x[(x$op=="=~"|x$op=="~")&is.na(x$ustart), ]$ustart <- 0
    }
    
    ## set the starting values in S and free parameters as 0 if NA
    if (any((x$op=="~~"&is.na(x$ustart)&x$free!=0))) {
        x[x$op=="~~"&is.na(x$ustart)&x$free!=0, ]$ustart <- 0
    }
        
    ## all variables
    ## Removed sort(); otherwise, the variables will be arranged as x1, x10, x2, x3...
    all.var <- unique(c(x$lhs, x$rhs))
    ## latent variables: (with indicators)
    latent <- unique(x[x$op== "=~", ]$lhs)
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
    Mmatrix <- matrix(0, nrow=1, ncol=no.all, dimnames=list(1, all.var))
  
    ## ## latent variable
    ## lhs <- (all.var %in% my.model[my.model$op== "=~", ]$lhs) 
    ## ## observed variables
    ## rhs <- (all.var %in% my.model[my.model$op== "=~", ]$rhs)

    ## Prepare the labels
    for (i in seq_len(nrow(x))) {
      ## if there is no label
      if (x[i, ]$label=="") {
        switch(x[i, ]$op,
               "=~" = x[i, ]$label <- paste0(x[i, ]$rhs, A.notation, x[i, ]$lhs),
               "~"  = x[i, ]$label <- paste0(x[i, ]$lhs, A.notation, x[i, ]$rhs),
               "~~" = x[i, ]$label <- paste0(x[i, ]$lhs, S.notation, x[i, ]$rhs),
               "~1" = x[i, ]$label <- paste0(x[i, ]$lhs, M.notation))
      }
    }
    
    ## replace NA to 0 in ustart
    x$ustart[is.na(x$ustart)] <- 0
    ## keys in as.mxMatrix format
    key <- with(x, ifelse(free==0, yes=ustart, no=paste(ustart, label, sep="*")))  
    
    for (i in seq_len(nrow(x))) {
        my.line <- x[i, ]
        switch(my.line$op,
               ## lhs: IV; rhs: DV
               "=~" = Amatrix[my.line$rhs, my.line$lhs] <- key[i],
               ## lhs: DV; rhs: IV
               "~" = Amatrix[my.line$lhs, my.line$rhs] <- key[i],
               "~~" = Smatrix[my.line$lhs, my.line$rhs] <- Smatrix[my.line$rhs, my.line$lhs] <- key[i],
               ## means
               "~1" = Mmatrix[1, my.line$lhs] <- key[i])
    }

    Fmatrix <- create.Fmatrix(c(rep(1, no.obs), rep(0, no.lat)), as.mxMatrix=FALSE)
    dimnames(Fmatrix) <- list(observed, all.var)

    list(A=Amatrix, S=Smatrix, F=Fmatrix, M=Mmatrix)
}


lavaan2RAM <- function(model, obs.variables = NULL, A.notation="ON", S.notation="WITH",
                       M.notation="mean", auto.var = TRUE, std.lv = TRUE,
                       ngroups=1, ...) {
    ## if (!requireNamespace("lavaan", quietly=TRUE))    
    ##     stop("\"lavaan\" package is required for this function.")

    ## Default: fix the latent independent variables at 1
    my.model <- lavaan::lavaanify(model, fixed.x = FALSE, auto.var=auto.var,
                                  std.lv=std.lv, ngroups=ngroups, ...)

    ## Drop rows with group==0 
    my.model <- my.model[!(my.model$group==0), ]
    
    if (max(my.model$group)==1) {
        out <- .lavaan2RAM(my.model, obs.variables=obs.variables, A.notation=A.notation,
                           S.notation=S.notation, M.notation=M.notation)
    } else {
        model.list <- split(my.model, my.model$group)
        out <- lapply(model.list, .lavaan2RAM,
                      obs.variables=obs.variables, A.notation=A.notation,
                      S.notation=S.notation, M.notation=M.notation)
    }    
    out
}
