
require(mnormt)
require(coda)
require(doParallel)
require(doRNG)
require(BDgraph)
require(foreach)


#devtools::install_github("mauroflorez/BMGM")
library(BMGM)

###### parallel
nrep <- 30
core_used <- 15
cl <- makeCluster(core_used)
registerDoParallel(cl)


################################################################################
#####------------------ Simulation Scenarios - To Change! -----------------#####

p = 10 #p = 20
n = 200 #1 cases: n = 200;
miss_rate  = 0 #2 cases: miss_rate = 0.1, miss_rate =0

################################################################################

out <- foreach(i=1:nrep, .options.RNG=123, .errorhandling = "pass") %dorng% {

  graph.sim <- BDgraph::bdgraph.sim(p = p, n = n, type = "mixed", graph = "random")

  true <- graph.sim$G
  X <- graph.sim$data

  type <- c(rep("d", 4), rep("c", 2), rep("m", 2), rep("c", 2))
  type2 <- type
  type2[type == "d"] <- "p"
  type2[type == "c"] <- "g"
  type2[type == "z"] <- "p"
  type2[type == "m"] <- "c"

  level2 <- rep(1, p)
  level2[type2== "c"] <- 2

  col_mis <- c()
  col_ctr <- c()
  #MAR:
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      if(true[i,j] == 1 && !(j %in% col_mis) && !(i %in% col_ctr)){
        col_mis <- append(col_mis, i)
        col_ctr <- append(col_ctr, j)
        break
      }
    }
  }

  ind_missing <- sample(1:length(col_mis), round(length(col_mis)/2))
  #Generating Missing Data (MAR)
  X_NA <- missMethods::delete_MAR_censoring(X, miss_rate, col_mis[ind_missing], col_ctr[ind_missing])
  #Complete Cases:
  X_comp <- na.omit(X_NA)

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

  fit.bdgraph <- BDgraph::bdgraph(X_NA, method = "gcgm",
                                  not.cont = (type2 != "g" )*1, cores=1)
  adj.bd <- fit.bdgraph$p_links
  signs.bd <- sign(fit.bdgraph$K_hat)

  prop.bdgraph <- adj.bd*signs.bd
  prop.bdgraph <- prop.bdgraph[upper.tri(prop.bdgraph)]

  est.bdgraph <- summary(fit.bdgraph)$selected_g
  est.bdgraph <- est.bdgraph[upper.tri(est.bdgraph)]

  #MGM
  fit.mgm <- mgm::mgm(X_comp, type2, level2, threshold = "none")
  adj.mgm <- fit.mgm$pairwise$wadj
  sgn.mgm <- fit.mgm$pairwise$signs
  sgn.mgm[is.na(sgn.mgm)] <- 0
  prop.mgm <- adj.mgm*sgn.mgm
  prop.mgm <- prop.mgm[upper.tri(prop.mgm)]

  est.mgm <- (mgm::mgm(X_comp, type2, level2, threshold = "LW")$pairwise$wadj > 0)*1
  est.mgm <- est.mgm[upper.tri(est.mgm)]

  #OURS

  fit.our <- BMGM::bmgm(X = X_NA, type = type, nburn = 10000, nsample = 40000, v_0 = 0.01, v_1 = 1, bfdr = 0.05)
  adj.fit <- colMeans(fit.our$post_G)
  sgn.fit <- -sign(colMeans(fit.our$post_Beta))

  prop.our <- c(adj.fit*sgn.fit)
  names(prop.our) <- NULL

  est.our <- (fit.our$adj_G)*1
  est.our <- est.our[upper.tri(est.our)]


  true <- true[upper.tri(true)]

  summary <- c(round(c(prop.huge, prop.bdgraph, prop.mgm, prop.our, true), 2), c(est.huge, est.bdgraph, est.mgm, est.bdgraph, est.our))
}

mat_ris <- do.call(rbind, out)
nomefile <- file.path(paste0("Simulation_bdgraph_p_",p,"_n_",n,"_mr_", miss_rate,".csv"))
write.csv(mat_ris, nomefile)
stopCluster(cl)
