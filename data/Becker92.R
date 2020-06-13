r <- matrix(c(.46,.31,.19,.46,.55,.32,.397,.402,.183,.266,.567,.218,.28,.19,.18,.47,-.21,-.15), ncol=3, byrow=T)
my.df <- lapply(split(r, 1:6), function(x) { out <- diag(rep(1,3))
                                             out[lower.tri(out)] <- x
                                             out[upper.tri(out)] <- x
                                             dimnames(out) <- list(c("Math", "Spatial", "Verbal"),
                                                                   c("Math", "Spatial", "Verbal"))
                                             out})
names(my.df) <- c("Berry (1957)", "Rosenberg (1981)", "Weiner 1 (1984)", "Weiner 2 (1984)", "Becker 1 (1978)", "Becker 2 (1978)")
Becker92 <- list(data=my.df, n=c(103,69,69,70,153,74))   
rm(r, my.df)
