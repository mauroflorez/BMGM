
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

N = 200
miss_rate = 0.2

################################################################################
#####------------------ Simulation Scenarios - To Change! -----------------#####
library(mgm)

p <- 10
type <- c("p", "p", "p", "p", "g", "g", "c", "c", "g", "g")
level <- c(1,1,1,1,1,1,2,2,1,1)

# Running Models:

out <- foreach(i=1:nrep, .options.RNG=123, .errorhandling = "pass") %dorng% {

  simulated_data <- BDgraph::bdgraph.sim(p = 10, type = "mixed", n = N)

  true_adj <- simulated_data$G

  data <- simulated_data$data


  col_mis <- c()
  col_ctr <- c()
  #MAR:

  for(i in 1:(p-1)){
    for(j in (i+1):p){
      if(true_adj[i,j] == 1 && !(j %in% col_mis) && !(i %in% col_ctr)){
        col_mis <- append(col_mis, i)
        col_ctr <- append(col_ctr, j)
        break
      }
    }
  }

  ind_missing <- sample(1:length(col_mis), round(length(col_mis)/2))
  #Generating Missing Data (MAR)
  X_NA <- missMethods::delete_MAR_censoring(data, miss_rate, col_mis[ind_missing], col_ctr[ind_missing])
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

  est.mgm <- (mgm::mgm(data = data, type = type, level = level, threshold = "LW")$pairwise$wadj > 0)*1
  est.mgm <- est.mgm[upper.tri(est.mgm)]

  #OURS

  type_bmgm <- c(rep("z", 4), rep("c", 2), rep("m", 2), "c", "c") # "g", "p", "c"
  # Fit BMGM explicitly
  fit.our <- bmgm(X = X_NA, type = type_bmgm, nburn = 5000, nsample = 10000, v_0 = 0.05, v_1 = 1, bfdr = 0.05)
  adj.fit <- colMeans(fit.our$post_G)
  sgn.fit <- -sign(colMeans(fit.our$post_Beta))

  prop.our <- c(adj.fit*sgn.fit)
  names(prop.our) <- NULL

  est.our <- (fit.our$adj_G)*1
  est.our <- est.our[upper.tri(est.our)]

  fit.our2 <- bmgm(X_NA, type = type_bmgm, nburn = 5000, nsample = 10000, v_0 = 0.1, v_1 = 2)
  adj.fit2 <- colMeans(fit.our2$post_G)
  sgn.fit2 <- -sign(colMeans(fit.our2$post_Beta))

  prop.our2 <- c(adj.fit2*sgn.fit2)
  names(prop.our2) <- NULL

  est.our2 <- (fit.our2$adj_G)*1
  est.our2 <- est.our2[upper.tri(est.our2)]

  true <- true[upper.tri(true)]

  summary <- c(round(c(prop.huge, prop.bdgraph, prop.mgm, prop.our, prop.our2, true), 2), c(est.huge, est.bdgraph, est.mgm, est.bdgraph, est.our, est.our2))

  summary
}

mat_ris <- do.call(rbind, out)
nomefile <- file.path(paste0("Simulation_bdgraph_p_",p,"_N_",N,"_mr_", miss_rate,".csv"))
write.csv(mat_ris, nomefile)
stopCluster(cl)

