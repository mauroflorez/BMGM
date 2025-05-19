#########################
##  Simulazioni Mauro
#########################

# R CMD BATCH prova_M.R out_M.txt &
# top -u d098100
# htop

#Require::setLibPaths("/home/d098100/R/x86_64-pc-linux-gnu-library/4.3", standAlone = FALSE)
#.libPaths()

require(mnormt)
require(coda)
require(doParallel)
require(doRNG)
require(foreach)
require(missMethods)

source("bmgm.R") #Load the main function
source("sampler_bmgm.R")
###### Parallel Cores
nrep <- 20
core_used <- 20
cl <- makeCluster(core_used)
registerDoParallel(cl)

################################################################################
#####------------------ Simulation Scenarios - To Change! -----------------#####

p = 20 #p = 20
n = 200 #1 cases: n = 200;
miss_rate  = 0 #3 cases: miss_rate = 0.1, miss_rate = 0.2

################################################################################

#Fixed Parameters and comparisons:
k <- p/2
type <- rep(c("c", "d"), k)
categories <- rep(1, p)
type2 <- rep(c("g", "p"), k)
level2 <- rep(1, p)
lambda <- 5
prob_p <- 2/(p-1)

# Additional Functions:

F_transformation <- function(X, type, parameter, cont = FALSE){
  X_t <- X
  transform <- function(x, typ, param){
    switch(typ,
           m = x,
           d = atan(x)^param,
           z = atan(x)^param,
           c = if(cont) sign(x)*atan(abs(x))^param else x)
  }

  if(!(is.matrix(X) | is.data.frame(X))) X_t <- transform(X, type, parameter) else
    for(c in 1:ncol(X)) X_t[,c] <- transform(X[,c], type[c], parameter)

  return(X_t)
}

# Running Models:

out <- foreach(i=1:nrep, .options.RNG=123, .errorhandling = "pass") %dorng% {
    # gen <- generate_multi(n = n, k = k, seed = i)

    theta <- list()
    for(i in 1:p){
      switch(type[i],
             m = {
               prop <- runif(1, min = 0.4, 0.6)
               tta <- c(prop, 1 - prop)
             },
             c = {
               sd <- sample(c(1, 2, 4), 1)
               mean <- 0
               tta <- c(mean, sd)
             },
             d = {
               mu <- 1
               nu <- sample(c(0.5, 1, 2), 1)
               tta <- c(mu, nu)
             },
             z = {
               prop <- runif(1, min = 0.4, 0.8)
               mu <- sample(c(2, 4, 6), 1)
               tta <- c(prop, mu)
             })
      theta[[i]] <- tta
    }

    Beta <- diag(p)

    pairs <- combn(p, 2, simplify = F)

    for(e in pairs){
      Beta[e[1], e[2]] = Beta[e[2], e[1]] = sample(c(0,1), 1, prob = c(1 - prob_p, prob_p))*sample(c(-1,1), 1)*runif(1, 0.2, 0.6)
    }

    X <- sampler_bmgm(n = n, Beta = Beta, theta = theta, type = type, categories = categories, lambda = 5, M = 1000)

    true <- (Beta != 0)*1

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

    #huge

    fit.huge <- huge::huge(X_comp, method = "glasso")
    matrices <- fit.huge$path
    sgn.huge <- sign(fit.huge$icov[[length(fit.huge$icov)]])
    adj.huge <- apply(array(unlist(matrices), dim = c(dim(matrices[[1]]), length(matrices))),
                      MARGIN = c(1, 2), FUN = mean)

    prop.huge <- adj.huge*sgn.huge

    prop.huge <- prop.huge[upper.tri(prop.huge)]

    #BDgraph

    fit.bdgraph <- BDgraph::bdgraph(X_NA, method = "gcgm",
                                    not.cont = (type2 != "g" )*1, cores=1)
    adj.bd <- fit.bdgraph$p_links
    signs.bd <- sign(fit.bdgraph$K_hat)

    prop.bdgraph <- adj.bd*signs.bd

    prop.bdgraph <- prop.bdgraph[upper.tri(prop.bdgraph)]

    #MGM
    fit.mgm <- mgm::mgm(X_comp, type2, level2, threshold = "none")
    adj.mgm <- fit.mgm$pairwise$wadj
    sgn.mgm <- fit.mgm$pairwise$signs
    sgn.mgm[is.na(sgn.mgm)] <- 0

    prop.mgm <- adj.mgm*sgn.mgm

    prop.mgm <- prop.mgm[upper.tri(prop.mgm)]

    #OURS

    fit.our <- bmgm(X = X_NA, type = type, nburn = 10000, nsample = 10000, v_0 = 0.05, v_1 = 1, bfdr = 0.05)
    adj.fit <- colMeans(fit.our$post_G)
    sgn.fit <- sign(colMeans(fit.our$post_Beta))

    prop.our <- c(adj.fit*sgn.fit)
    names(prop.our) <- NULL

    true <- true[col(true) > row(true)]

    summary <- round(c(prop.huge, prop.bdgraph, prop.mgm, prop.our, true), 2)

    summary
  }

mat_ris <- do.call(rbind, out)
nomefile <- file.path(paste0("Simulation_imputation_p_",p,"_n_",n,"_mr_", miss_rate,".csv"))
write.csv(mat_ris, nomefile)
stopCluster(cl)

