#' Bayesian Mixed Graphical Model
#' @importFrom stats cov dbeta dgamma dnorm dpois rbinom rnorm sd var
#' @description Implements a Bayesian approach for inference on mixed graphical models.
#' Estimates conditional independencies between variables of mixed types
#' (continuous, discrete, categorical, zero-inflated count data) using an MCMC
#' algorithm for posterior inference.
#'
#' @param X A numeric matrix or data frame containing the data (mixed type variables).
#' @param type A character vector indicating the type of each variable.
#'        Options: "c" (continuous), "d" (discrete), "z" (zero-inflated count), "m" (categorical).
#' @param nburn Number of burn-in iterations for MCMC. Default is 1000.
#' @param nsample Number of MCMC samples to collect after burn-in. Default is 1000.
#' @param theta_priors List of priors for node parameters.
#' @param v_0 Prior variance for the spike in spike-and-slab. Default is 0.05.
#' @param v_1 Prior variance for the slab in spike-and-slab. Default is 1.
#' @param pi_beta Prior probability of edge inclusion. Default is `2/(p-1)`.
#' @param seed Random seed for reproducibility.
#' @param context_spec Logical, whether to compute the context-specific adjacency matrix. Default is TRUE.
#' @param bfdr Bayesian False Discovery Rate threshold for edge selection. Default is 0.05.
#' @param cont Logical, whether to transform continuous variables. Default is FALSE.
#' @param ... Additional arguments.
#'
#' @return A list containing:
#' \describe{
#'   \item{post_Beta}{Posterior samples of the precision matrix.}
#'   \item{post_theta}{Posterior samples of node parameters.}
#'   \item{post_Z}{Posterior samples of edge indicators.}
#'   \item{adj_Beta}{Estimated adjacency matrix for dependencies.}
#'   \item{adj_G}{Adjacency matrix for edge selection.}
#'   \item{lambda}{Estimated transformation parameter.}
#' }
#'
#' @export
#'
#' @examples
#' # Simulated example data
#' set.seed(123)
#' n <- 100
#' X1 <- rpois(n, lambda = 3)
#' X2 <- rnorm(n, mean = 5)
#' X3 <- rnorm(n, mean = X2)
#' X4 <- rnorm(n, mean = X1)
#' X <- cbind(X1,X2,X3,X4)
#' type <- c("d", "c", "c", "c")
#'
#' # Fit Bayesian Mixed Graphical Model
#' fit <- bmgm(X, type, nburn = 1000, nsample = 1000)
#'
#' # Print adjacency matrix
#' print(fit$adj_G)
#'
bmgm <- function(X, type, nburn = 1000, nsample = 1000, theta_priors,
                   v_0 = 0.05, v_1 = 1, pi_beta, seed, context_spec = T,
                   bfdr = 0.05, cont = FALSE,...){

  if(!missing(seed)) set.seed(seed)

  if(missing(type)){
    stop("Error: type of variables is missing")
  }

  X_input <- X

  n <- nrow(X)
  #Number of predictors: p
  p <- ncol(X)

  #Initial Imputation
  R <- 1 - is.na(X)*1
  r_imp <- which(rowSums(R) < p)
  means <- colMeans(X, na.rm = T)

  if(length(r_imp) > 0) {
    for(s in 1:p){
      impute_value <- switch(type[s],
                             "c" = means[s],
                             "d" = round(means[s]),
                             "z" = round(means[s]),
                             "m" = as.numeric(names(which.max(table(X[, s])))))
      X[, s][R[, s] == 0] <- impute_value
    }
  }
      #switch(type[s],
      #       c = {
      #         mean = means[s]
      #       },
      #       d = {
      #         mean = round(means[s])
      #       },
      #      z = {
      #         mean = round(means[s])
      #       },
      #       m = {
      #         mean = as.numeric(levels(factor(X[,s]))[summary(factor(X[,s])) == max(summary(factor(X[,s])))])
  #       })
  #X[,s][R[,s] == 0] <- mean

  #Centering of continuous
  X[, type == "c"] <- scale(X[, type == "c"], center = means[type == "c"], scale = F)

  #Spike and slab
  log_prior_beta <- function(Beta, pi_beta, v0, v1){
    sum(log(pi_beta*exp(dnorm(Beta, mean = 0, sd = v1, log = TRUE)) +
              (1 - pi_beta)*exp(dnorm(Beta, mean = 0, sd = v0, log = TRUE))))
  }

  #Split categorical variables in columns (if any):
  split <- split_X_cat(X, type)
  X_design <- split$matrix
  categories <- split$categories

  #Number of nodes: q
  q <- ncol(X_design)
  #Type of each node
  type_q <- rep(type, categories)

  #############======= Auxiliary functions to get the names =======#############

  #Get the names of the nodes
  var_names <- rep(1:p, categories)
  tags <- get_names_graph(p, q, categories, var_names)

  #Get the subnames for categorical nodes
  tag <- tags$tag_nodes
  #Get the names for the edges
  tag_Beta <- tags$tag_Beta
  #Get the indicators where there are no edges (among sub-categories)
  ind_noedge <- tags$indicators_noedge

  ####################======= Transformation ========###########################

  lambda <- find_lambda(X, type)

  F_X <- F_transformation(X = X_design, type = type_q, parameter = lambda, cont)
  std_err <- apply(F_X, 2, sd)
  F_scaled <- scale(F_X, center = F, scale = std_err)

  #Beta
  F_centered <- scale(F_scaled, center = T, scale = F)
  S <- crossprod(F_centered)
  pdxid <- diag(solve(var(F_scaled)))

  ########################===== Hyperpriors =====##############################

  #Theta_priors should be a list where each element corresponds to the priors
  #for each variable.
  #In the case of categorical data simply the uniform. By the default are given
  #as below, should be specified in specific scenarios.

  theta_priors <- list()
  for (s in 1:p) {
    theta_priors[[s]] <- switch(type[s],
                                "c" = c(mean_mu = 0.001, mean_sd = 0.001, sd_shape = 0.001, sd_rate = 0.001),
                                "d" = c(mu_shape = 0.001, mu_rate = 0.001, nu_shape = 0.001, nu_rate = 0.001),
                                "z" = c(p_shape1 = 0.001, p_shape2 = 0.001, mu_shape = 0.001, mu_rate = 0.001),
                                "m" = rep(1 / (categories[s] + 1), categories[s] + 1))
  }

  ##### Priors for Beta / Now as inputs

  if(missing(pi_beta)) pi_beta <- 2/(p-1)

  Beta <- diag(pdxid)
  G <- diag(q) - diag(q)
  theta <- list()

  for(s in 1:p){
    theta[[s]] <- switch(type[s],
                         'c' = c(0,1),
                         'd' = c(0.1, 1),
                         'z' = c(0.5, 1),
                         'm' = rep(1/(categories[s]+1), categories[s]+1))
  }

  # Initial variance for steps in MCMC
  h_beta <- rep(0.1, q)
  h_theta <- list()
  scale_theta <- c()
  for(s in 1:p){
    h_theta[[s]] <- switch(type[s],
                           'c' = diag(0.1, 2),
                           'd' = diag(c(0.01, 0.01)),
                           'm' = diag(0.1, categories[s]+1),
                           'z' = diag(c(0.005, 0.5)))

    scale_theta[s] <- switch(type[s],
                             'c' = 2.4^2/2,
                             'd' = 2.4^2/2,
                             'z' = 2.4^2/2,
                             'm' = 2.4^2/(categories[s]+1))
  }

  ac_anterior = ac_Beta = ac_gamma = ac_theta = rep(0, q)

  ############################===== Storage =====###############################

  post_Beta <- matrix(nrow = nburn+nsample, ncol = q*(q-1)/2)
  post_G <- matrix(nrow = nburn+nsample, ncol = q*(q-1)/2)
  post_imputation <- list()
  post_theta <- list()

  for(s in 1:p){
    post_theta[[s]] <- switch(type[s],
                              'c' = matrix(nrow = nburn+nsample, ncol = 2),
                              'd' = matrix(nrow = nburn+nsample, ncol = 2),
                              'z' = matrix(nrow = nburn+nsample, ncol = 2),
                              'm' = matrix(nrow = nburn+nsample, ncol = categories[s]+1))
  }

  colnames(post_G) <- tag_Beta[upper.tri(tag_Beta)]
  colnames(post_Beta) <- tag_Beta[upper.tri(tag_Beta)]

  log_norm_constant <- function(s, parameters){
    domain <- 0:100
    se <- std_err[s]
    un_llk_z <- function(x, theta)
      exp(theta[2]*(log(theta[1])*x - lfactorial(x)) - theta[3]*F_transformation(x, type = "d", lambda)/se)

    log(apply(parameters, 1, function(x) sum(un_llk_z(domain, theta = x))))
  }

  log_norm_constant_Z <- function(s, parameters){
    domain <- 0:100
    se <- std_err[s]
    un_llk <- function(x, theta){
      (theta[1]*(x == 0) + (1 - theta[1])*dpois(x, lambda = theta[2]))*
        exp(-theta[3]*F_transformation(x, type = "d", lambda)/se)
    }
    log(apply(parameters, 1, function(x) sum(un_llk(domain, theta = x))))
  }

  pb <- progress::progress_bar$new(format = "Completing [:bar] :percent || ETA: :eta]",
                                   total = nburn + nsample,
                                   clear = FALSE)

  for(m in 1:(nburn + nsample)){
    ######################## 1st Block - Updates Theta ###########################
    for(s in 1:p){
      theta_s <- theta[[s]]
      h_theta_s <- h_theta[[s]]

      #Proposal
      switch(type[s],
             'c' = {
               #sampling mu
               mu0 = theta_priors[[s]][1]
               tau0 = theta_priors[[s]][2]

               mu = means[s]#mean(X[,s])
               tao = theta_s[2]

               var_post <- tau0 + n*tao
               mean_post <- (tau0*mu0 + n*tao*mu)/var_post

               #adding norm. constant
               C_s = sum(c(F_scaled[,-s]%*%Beta[s,-s]))

               mean_post <- mean_post - (1/var_post)*C_s
               theta_s[1] <- rnorm(1, mean = mean_post, sd = 1/sqrt(var_post))

               #sampling tau. In this case we need MCMC

               tau <- theta_s[2]
               mu <- theta_s[1]
               tau_proposal <- abs(rnorm(1, mean = tau, h_theta_s))

               a0 = theta_priors[[s]][3]
               b0 = theta_priors[[s]][4]

               shape_post = a0 + n/2
               rate_post = b0 + 0.5*sum((X[,s] - mu)^2)

               C_s_2 = sum(c(F_scaled[,-s]%*%Beta[s,-s])^2)

               ar <- dgamma(tau_proposal, shape = shape_post,rate = rate_post, log = T) -
                 dgamma(tau, shape = shape_post, rate = rate_post, log = T) +
                 C_s_2/2*(1/tau - 1/tau_proposal)

               accept <- min(1, exp(ar))

               if(stats::runif(1) < accept){
                 theta_s[2] <- tau_proposal
                 ac_theta[s] <- ac_theta[s] + 1
               }

             },
             'd' = {
               theta_star <- abs(MASS::mvrnorm(n = 1, mu = c(theta_s[1], theta_s[2]),
                                               Sigma = h_theta_s))

               C_s = c(F_scaled[,-s]%*%Beta[s,-s])

               param_s <- cbind('mu' = rep(theta_s[1], n), 'nu' = rep(theta_s[2], n), 'edge' = C_s)
               param_star <- cbind('mu' = rep(theta_star[1], n), 'nu' = rep(theta_star[2], n), 'edge' = C_s)

               log_Z_s <- sum(log_norm_constant(s, param_s))
               log_Z_star <- sum(log_norm_constant(s, param_star))

               log_density_d <- sum(theta_s[2]*(log(theta_s[1])*X[,s] - lfactorial(X[,s])))
               log_density_d_star <- sum(theta_star[2]*(log(theta_star[1])*X[,s] - lfactorial(X[,s])))

               theta_priors_s <- theta_priors[[s]]

               log_prior_d <- dgamma(theta_s[1], shape = theta_priors_s[1], rate = theta_priors_s[2], log = TRUE) +
                 dgamma(theta_s[2], shape = theta_priors_s[3], rate = theta_priors_s[4], log = TRUE)

               log_prior_d_star <- dgamma(theta_star[1], shape = theta_priors_s[1], rate = theta_priors_s[2], log = TRUE) +
                 dgamma(theta_star[2], shape = theta_priors_s[3], rate = theta_priors_s[4], log = TRUE)

               #acceptance:

               log_ar <- log_density_d_star - log_density_d +
                 log_prior_d_star - log_prior_d +
                 log_Z_s - log_Z_star

               accept <- min(1, exp(log_ar))

               if(stats::runif(1) < accept){
                 theta_s <- theta_star
                 ac_theta[s] <- ac_theta[s] + 1
               }

             },
             'z' = {
               theta_star <- c(mnormt::rmtruncnorm(n = 1, mean = theta_s[1], varcov = h_theta_s[1,1], lower = 0, upper = 1),
                               abs(rnorm(n = 1, mean = theta_s[2], sd = h_theta_s[2,2])))

               C_s <- c(F_scaled[,-s]%*%Beta[s,-s])

               prior_values <- theta_priors[[s]]
               alpha0 <- prior_values[1] #assc with p
               beta0 <- prior_values[2]

               a0 <- prior_values[3] #with mu
               b0 <- prior_values[4]

               m_Z = sum(X[,s] == 0)

               post_s1 <- m_Z + alpha0
               post_s2 <- n - m_Z + beta0

               post_alpha <- a0 + sum(X[,s])
               post_beta <- n - m_Z + b0

               #Normalizing constant

               param_s <- cbind('pi' = rep(theta_s[1], n), 'mu' = rep(theta_s[2], n), 'edge' = C_s)
               param_star <- cbind('pi' = rep(theta_star[1], n), 'mu' = rep(theta_star[2], n), 'edge' = C_s)

               log_Z_s <- sum(log_norm_constant_Z(s, param_s))
               log_Z_star <- sum(log_norm_constant_Z(s, param_star))

               log_ar <- dbeta(theta_star[1], post_s1, post_s2, log = T) - dbeta(theta_s[1], post_s1, post_s2, log = T) +
                 dgamma(theta_star[2], post_alpha, post_beta, log = T) - dgamma(theta_s[2], post_alpha, post_beta, log = T) +
                 log_Z_s - log_Z_star

               accept <- min(1, exp(log_ar))

               if(stats::runif(1) < accept){
                 theta_s <- theta_star
                 ac_theta[s] <- ac_theta[s] + 1
               }
             },
             'm' = {
               theta_star <- abs(MASS::mvrnorm(n = 1, mu = theta_s, Sigma = h_theta_s))
               theta_star <- theta_star/sum(theta_star)

               cat = as.numeric(factor(X[,s]))
               #edge-potentials
               cols = which(var_names == s)
               se <- std_err[cols]
               C_s <- F_scaled[,-cols]%*%Beta[-cols,cols]

               un_llk <- cbind(rep(exp(log(theta_s[1])),n),
                               exp(apply(C_s, 1, function(x) log(theta_s[-1]) - x/se)))
               norm_consts <- rowSums(un_llk)

               un_llk_star <-  cbind(rep(exp(log(theta_star[1])),n),
                                     (exp(apply(C_s, 1, function(x) log(theta_star[-1]) - x/se))))
               norm_consts_star <- rowSums(un_llk_star)

               llk  <- un_llk[cbind(1:nrow(C_s), cat)]/norm_consts
               llk_star  <- un_llk_star[cbind(1:nrow(C_s), cat)]/norm_consts

               log_dif_llk <- sum(log(llk_star) - log(llk))

               prior_s <- theta_priors[[s]]
               #priors
               log_dif_priors <- log(gtools::ddirichlet(theta_star, prior_s)) -
                 log(gtools::ddirichlet(theta_s, prior_s))

               log_ar <- log_dif_llk + log_dif_priors

               accept <- min(1, exp(log_ar))

               if(stats::runif(1) < accept){
                 theta_s <- theta_star
                 ac_theta[s] <- ac_theta[s] + 1
               }
             })

      theta[[s]] <- theta_s
      post_theta[[s]][m,] <- theta[[s]]
    }

    ########################### 2nd Block - Update Beta ##########################
    # Arkaprava & Dunson (2020)

    #update = "beta"
    Beta_star = Beta
    for(l in 1:q){
      theta_s <- theta[[var_names[l]]]
      mean <- S[l,-l] #mean
      Omega <- Beta[-l, -l]
      diag(Omega) <- pmax(pdxid[-l], 1e-6)
      Omegai <- eigen(Omega)
      OmegatempiU <- t(Omegai$vectors)/sqrt(abs(Omegai$values))
      #OmegaiU <- Omegai$vectors
      #OmegaiD <- Omegai$values

      #Update column l
      Omega_inv <- crossprod(OmegatempiU)

      #Ci <- eigen((var(F_X[,var_names[l]]) + 1)*n*Omega_inv + diag(ifelse(G[l,-l] == 0, 1/v_0, 1/v_1)))
      Ci <- eigen((S[l,l] + 1)*Omega_inv + diag(ifelse(G[l,-l] == 0, 1/v_0, 1/v_1)))
      CiU <- t(Ci$vectors)/sqrt(abs(Ci$values))
      C_inv <- crossprod(CiU)
      #C_inv <- chol2inv(chol((S[l,l] + 1)*Omega_inv + diag(ifelse(G[l,-l] == 0, 1/v_0, 1/v_1))))

      #Proposal
      mean_proposal <- -C_inv%*%mean
      var_proposal <- C_inv

      Beta_proposal <- MASS::mvrnorm(1, mean_proposal, var_proposal)
      Beta_proposal[which(is.na(Beta_proposal))] <- 0

      #Adjust the update wrt the acceptance rate
      k_2 <- min(1, as.numeric(h_beta[l]/sqrt(crossprod(Beta_proposal - Beta[l,-l]))))
      Beta_proposal <- Beta[l,-l] + k_2*(Beta_proposal - Beta[l, -l])

      Beta_star[l, -l] <- Beta_proposal
      Beta_star[-l, l] <- Beta_proposal
      Beta_star <- Beta_star*ind_noedge

      C_s = c(F_scaled[,-l]%*%Beta[l,-l])
      C_star = c(F_scaled[,-l]%*%Beta_proposal)

      if(type[var_names[l]] == "c"){
        mu = theta_s[1]
        tau = theta_s[2]

        mean = mu - C_s/tau
        mean_star = mu - C_star/tau

        log_dif_llk <- sum(dnorm(X[,l], mean = mean_star, sd = sqrt(1/tau), log = T)) -
          sum(dnorm(X[,l], mean = mean, sd = sqrt(1/tau), log = T))
      } else {
        log_dif_llk <- sum(F_scaled[,l]*C_star) - sum(F_scaled[,l]*C_s)
      }

      log_dif_priors <- log_prior_beta(Beta_proposal, pi_beta, v_0, v_1) -
        log_prior_beta(Beta[l,-l], pi_beta, v_0, v_1)

      log_dif_prop <- mvtnorm::dmvnorm(Beta[-l, l], mean_proposal, var_proposal, log = T) -
        mvtnorm::dmvnorm(Beta_proposal, mean_proposal, var_proposal, log = T)

      switch(type[var_names[l]],
             'd' = {
               param_s <- cbind('mu' = rep(theta_s[1], n),
                                'nu' = rep(theta_s[2], n),
                                'edge' = C_s)
               param_star <- cbind('mu' = rep(theta_s[1], n),
                                   'nu' = rep(theta_s[2], n),
                                   'edge' = C_star)

               log_Z_s <- sum(log_norm_constant(s, param_s))
               log_Z_star <- sum(log_norm_constant(s, param_star))
               log_dif_norm <- log_Z_s - log_Z_star
             },
             'c' = {
               mu = theta_s[1]
               tau = theta_s[2]
               log_Z_s <- sum(mu*tau*C_s - C_s^2/(2*tau))
               log_Z_star <- sum(mu*tau*C_star - C_star^2/(2*tau))
               log_dif_norm = log_Z_s - log_Z_star #<- sum(mu*tau*(C_s - C_star) + 1/(2*tau)*(C_s^2 - C_star^2))
             },
             'z' = {
               param_s <- cbind('p' = rep(theta_s[1], n),
                                'mu' = rep(theta_s[2], n),
                                'edge' = C_s)
               param_star <- cbind('p' = rep(theta_s[1], n),
                                   'mu' = rep(theta_s[2], n),
                                   'edge' = C_star)

               log_Z_s <- sum(log_norm_constant_Z(s, param_s))
               log_Z_star <- sum(log_norm_constant_Z(s, param_star))
               log_dif_norm <- log_Z_s - log_Z_star
             },
             'm' = {
               Beta_star <- Beta
               Beta_star[l, -l] = Beta_star[-l,l] = Beta_proposal
               cat = as.numeric(factor(X[,s]))
               #edge-potentials
               cols = which(var_names == s)
               se <- std_err[cols]
               C_s <- F_scaled[,-cols]%*%Beta[-cols,cols]
               C_star <- F_scaled[,-cols]%*%Beta_star[-cols,cols]

               un_llk <- cbind(rep(exp(log(theta_s[1])),n),
                               exp(apply(C_s, 1, function(x) log(theta_s[-1]) - x/se)))
               log_Z_s <- sum(log(rowSums(un_llk)))

               un_llk_star <-  cbind(rep(exp(log(theta_s[1])),n),
                                     exp(apply(C_star, 1, function(x) log(theta_s[-1]) - x/se)))
               log_Z_star <- sum(log(rowSums(un_llk_star)))
               log_dif_norm <- log_Z_s - log_Z_star
             })

      log_ar <- log_dif_llk + log_dif_priors + log_dif_norm #+ log_dif_prop

      accept <- min(1, exp(log_ar))

      if(stats::runif(1) < accept){
        Beta = Beta_star
        ac_Beta[l] <- ac_Beta[l] + 1
      }
    }

    ######################== 3rd Block - Update G and pi ===###################
    G_new <- Beta[upper.tri(Beta)]
    slab <- pi_beta*dnorm(G_new, 0, v_1)
    spike <- (1-pi_beta)*dnorm(G_new, 0, v_0)
    G_new <- slab/(slab+spike)
    nan_ind <- is.nan(G_new)
    G_new[nan_ind] <- 0

    G_vector <- rbinom(q*(q-1)/2, size = 1, prob = G_new)

    G <- matrix(0, ncol = q, nrow = q)
    G[upper.tri(Beta)] <- G_vector
    G <- G + t(G)

    post_Beta[m,] <- Beta[upper.tri(Beta)]
    post_G[m,] <- G_vector


    #######################=== Adjust rates ===##############################

    if((m%%100)==0){# & m < nburn){
      ar_beta <- ac_Beta / m

      ar_theta <- ac_theta / m

      #ar_beta_m <- (ac_Beta - ac_anterior)/100

      h_beta[ar_beta < .2] <- h_beta[ar_beta < .2]/2
      h_beta[ar_beta > .6] <- h_beta[ar_beta > .6]*2

      scale_theta[ar_theta < .2] <- scale_theta[ar_theta < .2]/2
      scale_theta[ar_theta > .6] <- scale_theta[ar_theta > .6]*2

      for(l in 1:p){
        h_theta[[l]] <- scale_theta[l]*cov(post_theta[[l]][1:m,]) +
          scale_theta[l]*0.00001*diag(length(theta[[l]]))
      }

      #Imputation
      for(i in r_imp){
        new_values <- sampler_bmgm(n = 10, Beta = Beta, theta = theta, type = type,
                                   categories = categories, lambda = lambda, std = std_err,
                                   X_new = X[i,], variables = c(1 - R[i,]))
        new_values_mean <- colMeans(new_values)

        # Default: use mean for all
        X[i,] <- new_values_mean

        # Now fix "m" type (categorical)
        for (j in 1:p) {
          if (type[j] == "m" && R[i, j] == 0) {
            # Take the mode (most common category)
            X[i, j] <- as.numeric(names(which.max(table(new_values[, j]))))
          }
        }

        # Then apply rounding for discrete, zero-inflated
        X[i, type %in% c("d", "z") & R[i,] == 0] <- round(X[i, type %in% c("d", "z") & R[i,] == 0])

        # Finally recenter continuous variables
        # not longer necessary as centering is done prior
        #X[i, type == "c" & R[i,] == 0] <- X[i, type == "c" & R[i,] == 0] - means[type == "c" & R[i,] == 0]
      }

      post_imputation[[m/100]] <- X

      split <- split_X_cat(X, type)
      X_design <- split$matrix
      lambda <- find_lambda(X, type)

      F_X <- F_transformation(X = X_design, type = type_q, parameter = lambda, cont)

      std_err <- apply(F_X, 2, sd)
      F_scaled <- scale(F_X, center = F, scale = std_err)

      #S <- crossprod(F_centered)
      #pdxid <- diag(solve(var(F_scaled)))
    }
    pb$tick()
  }

  ce_graph <- context_spec_graph(q, post_Beta, post_G, tag, bfdr)
  #General Graph
  cat_graph <- categories_graph(q, p, var_names, ce_graph$ce_esti_Z,
                                ce_graph$ce_esti_Beta, categories)

  fit <- list(post_Beta = post_Beta, post_theta = post_theta, post_G = post_G,
              adj_Beta = cat_graph$Adj_Beta, adj_G = (cat_graph$Adj_Z!=0)*1,
              lambda = lambda, std = std_err, X = X_input, type = type)

  if(context_spec == T & any(type == "m")){
    fit[["adj_Beta_ce"]] <- ce_graph$ce_esti_Beta
    fit[["adj_Z_ce"]] <- ce_graph$ce_esti_Z
  }

  return(fit)
}
