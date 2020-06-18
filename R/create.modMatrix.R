create.modMatrix <- function(RAM, output=c("A", "S"), mod) {
    output <- match.arg(output)

    switch(output,
        A = { out <- RAM$A
        out[grep("\\*", out)] <- paste0("0*data.", mod)},
        S = { out <- RAM$S
        out[grep("\\*", out)] <- paste0("0*data.", mod)
        Diag(out) <- "0"})

    out
}

