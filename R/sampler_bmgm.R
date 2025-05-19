#' Gibbs Sampler for Bayesian Mixed Graphical Models
#'
#' @description Implements a Gibbs sampler for Bayesian inference on Mixed Graphical Models.
#'
#' @param n Number of samples to generate.
#' @param Beta Precision matrix of the graphical model.
#' @param theta A list of parameters for each node.
#' @param type Character vector indicating variable types.
#' @param categories Vector indicating the number of categories (for categorical variables).
#' @param lambda Transformation parameter.
#' @param M Number of Gibbs sampler iterations. Default is 1000.
#' @param X_new Vector of known values for imputation. Default is `c(0,0,0)`.
#' @param variables Indicator vector specifying which variables to update.
#' @param std Standard deviation used in the transformation.
#'
#' @return A matrix containing the sampled data.
#'
#' @examples
#' # Simulate example data
#' set.seed(123)
#' Beta <- diag(4)  # Identity precision matrix
#' theta <- list(c(0,1), c(0,1), c(0.5,1), c(0.5,1))
#' type <- c("c", "c", "d", "d")
#' categories <- c(1,1,1,1)
#' lambda <- 1
#'
#' # Run the Gibbs sampler
#' samples <- sampler_bmgm(n = 10, Beta = Beta, theta = theta, type = type,
#'                         categories = categories, lambda = lambda, M = 500)
#'
#' # View first few samples
#' head(samples)
#'
#' @export
sampler_bmgm <- function(n = 1, Beta, theta, type, categories, lambda,
                         M = 1000, X_new, variables, std){

  design <- function(x, type, categories){
    p <- length(x)
    x_design <- c()

    for(s in 1:p){
      dec <- decompose(x[s], type[s], categories[s])
      x_design <- append(x_design, dec)
    }
    x_design
  }
  decompose <- function(x, type, cat) {
    # use switch to evaluate the type and return a vector
    if(type == "m"){
      x_new <- rep(0, cat)
      x_new[x] <- 1
      x_new
    } else x
  }

  if(any(type == "m") && missing(categories)) stop("Error: categories vector is needed.")
  p <- length(theta)
  k <- ncol(Beta)

  x_new <- rep(0, p) #node values
  if(missing(X_new)){
    X_new <- rep(0, k) #design matrix
  } else if(length(X_new) == p){
    x_new <- X_new
    X_new <- design(X_new, type, categories)
  }

  if(missing(variables)) variables <- rep(1, p)

  if(missing(categories)){
    categories = cat = rep(1, p)
  } else {
    cat = categories
    cat[type == "m"] <- categories[type == "m"] + 1
  }
  if(sum(categories) != k) stop("Dimensions are not matching!")

  type_k <- rep(type, categories)
  var_names <- rep(1:p, categories)

  #Matrix to save results
  x_sam <- matrix(nrow = M, ncol = k)
  x <- matrix(nrow = M, ncol = p)
  x_sam[1,] <- X_new
  x[1,] <- x_new


  #vector of updateable variables
  update <- which(variables == 1)
  possible <- (1:k)[var_names %in% update]
  #Scaling given or not
  ind = 0
  if(missing(std)){
    std <- rep(1, k)
    ind = 1
  }

  #Gibbs Sampler
  for(m in 2:M){
    x_sam[m, (1:k)[-possible]] <- x_sam[m-1, (1:k)[-possible]]
    x[m, (1:p)[-update]] <- x[m-1, (1:p)[-update]]

    for(i in update){
      l <- which(var_names == i)
      Beta_ts <- Beta[-l, l]
      tta <- theta[[i]]

      #Vector of updated values for interactions
      if(i == 1){
        aux <- c(x_sam[m-1, -l])
      } else if(i == p){
        aux <- c(x_sam[m, -c(min(l),k)])
      } else {
        aux <- c(x_sam[m, 1:(min(l)-1)], x_sam[m-1, (max(l)+1):k])
      }

      se <- std[l]

      edge_pot <- c((F_transformation(matrix(aux, nrow = 1),
                                      type_k[-i], lambda)/std[-l])%*%Beta_ts)

      switch(type[i],
             c = {
               #log_density_c <- function(x){
                #   dnorm(x, mean = tta[1], sd = tta[2], log = T) -
                #     edge_pot*(F_transformation(x, type = "c", parameter = lambda)/se)
               #}
               #x_star <- armspp::arms(1, log_density_c, -1000, 1000, metropolis = TRUE)
               #x[m,i] <- x_star
               mu_s <- tta[1]
               tau_s <- tta[2]

               mean_star <- mu_s - (1/tau_s) * (edge_pot / se)
               sd_star <- 1 / sqrt(tau_s)

               x_star <- rnorm(1, mean = mean_star, sd = sd_star)
               x[m,i] <- x_star
             },
             d = {
               log_density_d <- function(x){
                 tta[2]*(log(tta[1])*x - lfactorial(x)) -
                   edge_pot*(F_transformation(x, type = "d", parameter = lambda)/se)
               }

               x_star <- round(armspp::arms(1, log_density_d, 0, 1000, metropolis = TRUE))
               x[m,i] <- x_star
             },
             z = {
               log_density_z <- function(x){
                 log(tta[2])*x - lfactorial(x) -
                   edge_pot*(F_transformation(x, type = "z", parameter = lambda)/se)
               }

               x_star_s <- round(armspp::arms(1, log_density_z, 0, 1000, metropolis = TRUE))
               zero <- sample(c(0,1), 1, replace = T, prob = c(tta[1], 1 - tta[1]))

               x_star <- x_star_s*zero
               x[m,i] <- x_star
             },
             m = {
               #for each category i can calculate the likelihood
               ncat <- cat[i]
               x_star <- rep(0, ncat - 1)
               edge_pot <- append(0, edge_pot) #Insert category base
               se <- append(1, se)
               lik <- rep(0, ncat)
               for(c in 1:ncat){
                 lik[c] <- exp(log(tta[c]) - edge_pot[c]/se[c])
               }
               cat_sam <- sample(0:(ncat-1), 1, prob = lik/sum(lik))
               x[m,i] <- cat_sam
               x_star[cat_sam] <- 1
             })

      x_sam[m,l] <- x_star
    }
    if(ind == 1) std <- apply(F_transformation(x_sam[1:m,], type_k, lambda), 2, sd)
    std[is.na(std)] <- 1
    std[std == 0] <- 1

  }
  X_r <- x[sample(round(M/2):M, n),]
  return(X_r)
}
