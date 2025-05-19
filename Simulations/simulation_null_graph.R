
require(mnormt)
require(coda)
require(doParallel)
require(doRNG)
require(BDgraph)
require(foreach)

source("/home/mf53/bmgm.R")
source("/home/mf53/sampler_bmgm.R")
###### parallel
nrep <- 20
core_used <- 20
cl <- makeCluster(core_used)
registerDoParallel(cl)

N = 500
miss_rate = 0.2


################################################################################
#####------------------ Simulation Scenarios - To Change! -----------------#####
library(mgm)

simulate_null_graph_data <- function(n = 500, p = 10, type = c(rep("g",4), rep("p",3), rep("c",3)), level = c(1,1,1,1,1,1,1,2,2,2), missing_rate = 0.0) {
  X <- matrix(NA, n, p)

  for (j in 1:p) {
    if (type[j] == "g") {
      X[, j] <- rnorm(n, mean = 0, sd = 1)
    } else if (type[j] == "p") {
      X[, j] <- rpois(n, lambda = 4)
    } else if (type[j] == "c") {
      X[, j] <- sample(1:level[j], size = n, replace = TRUE) - 1
    }
  }

  # Add missingness if needed
  if (missing_rate > 0) {
    missing_mask <- matrix(rbinom(n*p, 1, missing_rate), n, p)
    X[missing_mask == 1] <- NA
  }

  colnames(X) <- paste0("V", 1:p)
  return(X)
}


p <- 10
type <- c("g", "g", "g", "g", "p", "p", "p", "c", "c", "c")
level <- c(1,1,1,1,1,1,1,2,2,2)


# Running Models:

out <- foreach(i=1:nrep, .options.RNG=123, .errorhandling = "pass") %dorng% {
  # gen <- generate_multi(n = n, k = k, seed = i)

  simulated_data <- simulate_null_graph_data(n = N, p = 10, missing_rate = miss_rate)
  true_adj <- matrix(0, p, p)

  colnames(true_adj) <- rownames(true_adj) <- paste0("V", 1:10)

  data <- simulated_data

  X_NA <- data
  #Complete Cases:
  X_comp <- na.omit(X_NA)

  not.cont = (type != "g")*1

  true <- true_adj
  #huge

  fit.huge <- huge::huge(X_comp, method = "glasso")
  matrices <- fit.huge$path
  sgn.huge <- sign(fit.huge$icov[[length(fit.huge$icov)]])
  adj.huge <- apply(array(unlist(matrices), dim = c(dim(matrices[[1]]), length(matrices))),
                    MARGIN = c(1, 2), FUN = mean)

  prop.huge <- adj.huge*sgn.huge
  prop.huge <- prop.huge[upper.tri(prop.huge)]

  est.huge <- huge::huge.select(fit.huge)$refit
  est.huge <- est.huge[upper.tri(est.huge)]

  #BDgraph
  not.cont = (type != "g")*1
  fit.bdgraph <- BDgraph::bdgraph(X_NA, method = "gcgm",
                                  not.cont = not.cont, cores=1)
  adj.bd <- fit.bdgraph$p_links
  signs.bd <- sign(fit.bdgraph$K_hat)

  prop.bdgraph <- adj.bd*signs.bd
  prop.bdgraph <- prop.bdgraph[upper.tri(prop.bdgraph)]

  est.bdgraph <- summary(fit.bdgraph)$selected_g
  est.bdgraph <- est.bdgraph[upper.tri(est.bdgraph)]

  #MGM
  fit.mgm <- mgm::mgm(data = X_comp, type = type, level = level, threshold = "none")
  adj.mgm <- fit.mgm$pairwise$wadj
  sgn.mgm <- fit.mgm$pairwise$signs
  sgn.mgm[is.na(sgn.mgm)] <- 0
  prop.mgm <- adj.mgm*sgn.mgm
  prop.mgm <- prop.mgm[upper.tri(prop.mgm)]

  est.mgm <- (mgm::mgm(data = X_comp, type = type, level = level, threshold = "LW")$pairwise$wadj > 0)*1
  est.mgm <- est.mgm[upper.tri(est.mgm)]

  #OURS

  type_bmgm <- c(rep("c", 4), rep("d", 3), rep("m", 3)) # "g", "p", "c"
  # Fit BMGM explicitly
  fit.our <- bmgm(X = X_NA, type = type_bmgm, nburn = 5000, nsample = 10000, v_0 = 0.05, v_1 = 1, bfdr = 0.05)
  adj.fit <- colMeans(fit.our$post_G)
  sgn.fit <- -sign(colMeans(fit.our$post_Beta))

  prop.our <- c(adj.fit*sgn.fit)
  names(prop.our) <- NULL

  est.our <- (fit.our$adj_G)*1
  est.our <- est.our[upper.tri(est.our)]

  fit.our2 <- bmgm(X_NA, type = type_bmgm, nburn = 10, nsample = 10, v_0 = 0.1, v_1 = 2)
  adj.fit2 <- colMeans(fit.our2$post_G)
  sgn.fit2 <- sign(colMeans(fit.our2$post_Beta))

  prop.our2 <- c(adj.fit2*sgn.fit2)
  names(prop.our2) <- NULL

  est.our2 <- (fit.our2$adj_G)*1
  est.our2 <- est.our2[upper.tri(est.our2)]

  true <- true[upper.tri(true)]

  summary <- c(round(c(prop.huge, prop.bdgraph, prop.mgm, prop.our, prop.our2, true), 2), c(est.huge, est.bdgraph, est.mgm, est.bdgraph, est.our, est.our2))

  summary
}

mat_ris <- do.call(rbind, out)
nomefile <- file.path(paste0("Simulation_null_p_",p,"_N_",N,"_mr_", miss_rate,".csv"))
write.csv(mat_ris, nomefile)
stopCluster(cl)

