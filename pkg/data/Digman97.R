Digman97 <- matrix(
c(-.48,-.10,.62,.27,.41,.59,.37,.00,.35,.41,
-.30,.07,.39,.09,.53,.59,.45,-.05,.44,.22,
.25,-.10,.65,.24,.35,.37,.41,.14,.33,.41,
-.26,-.16,.65,.01,.70,.71,.66,-.03,.24,.11,
.29,.16,.64,.32,.35,.27,.53,.22,.22,.36,
.35,.20,.66,.49,.57,.45,.59,.38,.31,.31,
.13,.43,.25,.37,.59,.28,.35,.15,.12,.10,
.16,.26,.36,.36,.41,.26,.33,.19,.16,.07,
.11,.19,.18,.22,.44,.42,.56,.24,.05,.12,
.42,.25,.34,.26,.69,.43,.46,.44,.54,.42,
.04,.27,.24,.21,.25,.53,.40,-.02,-.02,-.02,
-.07,.22,.13,.21,.25,.49,.43,-.06,-.04,-.05,
-.04,-.03,.25,-.03,.34,.41,.28,-.17,.08,.12,
.06,.04,.13,.16,.23,.17,.24,-.09,-.03,-.01),
ncol=10, byrow=TRUE)
Digman97 <- lapply(split(Digman97, 1:14), 
                 function(x) {mat <- matrix(1, ncol=5, nrow=5);
                              mat[upper.tri(mat, diag=FALSE)] <- x;
                              mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)];
                              mat})
Digman97 <- lapply(Digman97, function(x, var.names) {dimnames(x) <- list(var.names, var.names); x},
                 var.names=c("E", "A", "C", "ES", "I"))
names(Digman97) <- c("Digman 1 (1994)", "Digman 2 (1994)", "Digman 3 (1963c)", "Digman & Takemoto-Chock (1981b)",
           "Graziano & Ward (1992)", "Yik & Bond (1993)", "John et al. 1 (1984)", "John et al. 2 (1984)",
           "Costa & McCrae 1 (1992c)", "Costa & McCrae 2 (1992b)", "Costa & McCrae 3 (1992b)",
           "Costa, McCrae, & Dye (1991)", "Barrick & Mount (1993)", "Goldberg (1992a)")

Digman97.n <- c(102,149,334,162,91,656,70,70,277,227,1000,227,91,1040) 

Digman97.cluster <- c(rep("Children", 4), "Adolescents", rep("Young adults", 3), rep("Mature adults", 6))

