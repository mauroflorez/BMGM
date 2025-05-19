
require(mnormt)
require(coda)
require(doParallel)
require(doRNG)
require(foreach)
require(missMethods)

source("/home/mf53/bmgm.R")
source("/home/mf53/sampler_bmgm.R")

###### Parallel Cores
nrep <- 15
core_used <- 15
cl <- makeCluster(core_used)
registerDoParallel(cl)

################################################################################
#####------------------ Simulation Scenarios - To Change! -----------------#####

set.seed(123)
p = 12
n = 500
miss_rate  = 0 #0.1, 0.2

hyper_param <- matrix(c(0.01, 1, 2/(p-1),
                 0.025, 1, 2/(p-1),
                 0.05, 1, 2/(p-1),
                 0.075, 1, 2/(p-1),
                 0.1, 1, 2/(p-1),
                 0.05, 0.05, 2/(p-1),
                 0.05, 0.1, 2/(p-1),
                 0.05, 0.5, 2/(p-1),
                 0.05, 1, 2/(p-1),
                 0.05, 2, 2/(p-1),
                 0.05, 1, 1/(p-1),
                 0.05, 1, 2/(p-1),
                 0.05, 1, 3/(p-1),
                 0.05, 1, 4/(p-1),
                 0.05, 1, 5/(p-1)), ncol = 3, byrow = TRUE)


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

# Sample Data:

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

# Data:
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

out <- foreach(i=1:nrep, .options.RNG=123, .errorhandling = "pass") %dorng% {
  v_0 = hyper_param[i, 1]
  v_1 = hyper_param[i, 2]
  pi_beta = hyper_param[i, 3]
  #fitting
  t1<-Sys.time()
  fit.our <- bmgm(X = X_NA, type = type, nburn = 5000, nsample = 10000,
                  v_0 = v_0, v_1 = v_1, pi_beta = pi_beta)
  t2<-Sys.time()
  adj.fit <- (fit.our$adj_G>0)*1
  total_edges <- sum(adj.fit)/2
  comparison <- c(BDgraph::compare(adj.fit, true)[,2], total_edges)

  comparison
}

mat_ris <- do.call(rbind, out)
nomefile <- file.path(paste0("Simulation_sensitivity_p",p,"_n_",n,"_mr_", miss_rate,".csv"))
write.csv(mat_ris, nomefile)
stopCluster(cl)
