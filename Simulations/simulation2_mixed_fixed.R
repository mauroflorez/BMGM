
require(mnormt)
require(coda)
require(doParallel)
require(doRNG)
require(BDgraph)
require(foreach)

source("/home/mf53/bmgm.R")

###### parallel
nrep <- 30
core_used <- 15
cl <- makeCluster(core_used)
registerDoParallel(cl)


################################################################################
#####------------------ Simulation Scenarios - To Change! -----------------#####

p = 14 #p = 20
n = 500 #1 cases: n = 200;
miss_rate  = 0.1 #2 cases: miss_rate = 0.1, miss_rate =0

################################################################################

set.seed(123)

pairs <- list(
  c(1,2),
  c(2,3),
  c(2,4),
  c(4,5),
  c(3,6),
  c(7,8),
  c(9,10),
  c(9,12),
  c(4,12),
  c(3,10),
  c(5,7),
  c(5,9),
  c(14,5),
  c(13,7)
)

type <- c("m", "c", "z", "c", "z", "c", "z", "m", "c", "z", "c", "z", "c", "z")
level2 <- c(2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1)

type2 <- type
type2[type == "d"] <- "p"
type2[type == "c"] <- "g"
type2[type == "z"] <- "p"
type2[type == "m"] <- "c"

#Fixed Parameters and comparisons:
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
           mu <- 6
           nu <- sample(c(0.5, 1, 2), 1)
           tta <- c(mu, nu)
         },
         z = {
           prop <- runif(1, min = 0.4, 0.8)
           mu <- sample(c(1, 2, 5), 1)
           tta <- c(prop, mu)
         })
  theta[[i]] <- tta
}

Beta <- diag(p)

for(e in pairs){
  Beta[e[1], e[2]] = Beta[e[2], e[1]] = sample(c(-1,1), 1)*runif(1, 0.4, 0.8)
}

lambda <- 5

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

    X <- sampler_bmgm(n = n, Beta = Beta, theta = theta, type = type, categories = rep(1, p), lambda = 5, M = 2000)

    true <- (Beta != 0)*1

    col_mis <- sample(1:p, size = p/2)
    col_ctr <- sample((1:p)[!((1:p) %in% col_mis)], size = p/2)
    #MAR:
    #for(i in 1:(p-1)){
    #  for(j in (i+1):p){
    #    if(true[i,j] == 1 && !(j %in% col_mis) && !(i %in% col_ctr)){
    #      col_mis <- append(col_mis, i)
    #      col_ctr <- append(col_ctr, j)
    #      break
    #    }
    #  }
    #}

    #Generating Missing Data (MAR)
    X_NA <- missMethods::delete_MAR_censoring(X, miss_rate, col_mis, col_ctr)
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

    fit.our <- bmgm(X = X_NA, type = type, nburn = 5000, nsample = 5000, v_0 = 0.05, v_1 = 1, bfdr = 0.05)
    adj.fit <- colMeans(fit.our$post_G)
    sgn.fit <- sign(colMeans(fit.our$post_Beta))

    prop.our <- c(adj.fit*sgn.fit)
    names(prop.our) <- NULL

    est.our <- (fit.our$adj_Z)*1
    est.our <- est.our[upper.tri(est.our)]

    true <- true[upper.tri(true)]

    summary <- c(round(c(prop.huge, prop.bdgraph, prop.mgm, prop.our, true), 2), c(est.huge, est.bdgraph, est.mgm, est.bdgraph, est.our))

    summary
  }

mat_ris <- do.call(rbind, out)
nomefile <- file.path(paste0("Simulation2_fixed_imputation_p_",p,"_n_",n,"_mr_", miss_rate,".csv"))
write.csv(mat_ris, nomefile)


