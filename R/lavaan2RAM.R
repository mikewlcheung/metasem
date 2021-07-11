lavaan2RAM <- function(model, obs.variables = NULL, A.notation="ON", S.notation="WITH",
                       M.notation="mean", auto.var = TRUE, std.lv = TRUE,
                       ngroups=1, ...) {
    ## if (!requireNamespace("lavaan", quietly=TRUE))    
    ##     stop("\"lavaan\" package is required for this function.")

    ## Default: fix the latent independent variables at 1
    my.model <- lavaan::lavaanify(model, fixed.x = FALSE, auto.var=auto.var,
                                  std.lv=std.lv, ngroups=ngroups, ...)

    ## Maximum no. of groups
    max.gp <- max(my.model$group)

    ## Empty list to store the matrices per group
    out <- list()
    
    for (gp in seq_len(max.gp)) {
        mod <- my.model[my.model$group==gp, ]

        ## set the starting values in A and M as 0 if NA
        if (any((mod$op=="=~"|mod$op=="~"|mod$op=="~1")&is.na(mod$ustart))) {
            mod[(mod$op=="=~"|mod$op=="~"|mod$op=="~1")&is.na(mod$ustart), ]$ustart <- 0
        }
    
        ## set the starting values in S and free parameters as 0 if NA
        if (any((mod$op=="~~"&is.na(mod$ustart)&mod$free!=0))) {
            mod[mod$op=="~~"&is.na(mod$ustart)&mod$free!=0, ]$ustart <- 0
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
        Mmatrix <- matrix(0, nrow=1, ncol=no.all, dimnames=list(1, all.var))

        ## Prepare the labels
        for (i in seq_len(nrow(mod))) {
            ## if there is no label
            if (mod[i, ]$label=="") {
                switch(mod[i, ]$op,
                       "=~" = mod[i, ]$label <- paste0(mod[i, ]$rhs, A.notation, mod[i, ]$lhs),
                       "~"  = mod[i, ]$label <- paste0(mod[i, ]$lhs, A.notation, mod[i, ]$rhs),
                       "~~" = mod[i, ]$label <- paste0(mod[i, ]$lhs, S.notation, mod[i, ]$rhs),
                       "~1" = mod[i, ]$label <- paste0(mod[i, ]$lhs, M.notation))
            }
        }
    
        ## replace NA to 0 in ustart
        mod$ustart[is.na(mod$ustart)] <- 0
        ## keys in as.mxMatrix format
        key <- with(mod, ifelse(free==0, yes=ustart, no=paste(ustart, label, sep="*")))  
        
        for (i in seq_len(nrow(mod))) {
            my.line <- mod[i, ]
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
                                                   "<- mxAlgebra(", y[i,'rhs'], ", name=\"", y[i,'lhs'], "\")")))
                eval(parse(text=paste0("mxalgebra <- c(mxalgebra, ", y[i, 'lhs'], "=", y[i, 'lhs'], ")")))  },
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
        out[[1]]  <- list(A=out[[1]]$A, S=out[[1]]$S, F=out[[1]]$F, M=out[[1]]$M, mxalgebra=mxalgebra)
    }
    ## else {
    ##     out[[1]]  <- list(A=out[[1]]$A, S=out[[1]]$S, F=out[[1]]$F, M=out[[1]]$M)
    ## }


    ## Output the first list instead of a list of one item when there is only 1 group
    if (max.gp==1) {
        out  <- out[[1]]
    }
    
    out
}
