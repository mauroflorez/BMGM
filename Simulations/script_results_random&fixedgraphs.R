
library(corrplot)
library(ggcorrplot)
library(cowplot)
library(ggplot2)
library(latex2exp)

setwd("C:/Users/User/Documents")

p = 14
l = p*(p-1)/2
mr1 = 0
mr2 = 0.1
mr3 = 0.2

n1 = 500
n2 = 500
n3 = 500

#results1 <- read.csv(paste0("C:/Users/User/Documents/Simulation_imputation_p_",p,"_n_",n1,"_mr_",mr1,".csv"))[,-1]
results1 <- read.csv(paste0("C:/Users/User/Documents/Simulation2_fixed_imputation_new_p_",p,"_n_",n1,"_mr_",mr1,".csv"))[,-1]
results2 <- read.csv(paste0("C:/Users/User/Documents/Simulation2_fixed_imputation_new_p_",p,"_n_",n2,"_mr_",mr2,".csv"))[,-1]
results3 <- read.csv(paste0("C:/Users/User/Documents/Simulation2_fixed_imputation_new_p_",p,"_n_",n3,"_mr_",mr3,".csv"))[,-1]

results1 <- read.csv("C:/Users/User/Documents/Simulation_null_p_10_N_200_mr_0.1.csv")
results2 <- read.csv("C:/Users/User/Documents/Simulation_null_p_10_N_200_mr_0.2.csv")

#evaluation function at the end

#results1 <- as.data.frame(sim[,2:226])
results1 <- results1[,1:(l*5)]
results2 <- results2[,1:(l*5)]
results3 <- results3[,1:(l*5)]
ev1 <- evaluation(results1, n = n1, p = p, mr = mr1)
ev2 <- evaluation(results2, n = n2, p = p, mr = mr2)
ev3 <- evaluation(results3, n = n3, p = p, mr = mr3)

results_all <- rbind(ev1$results,
                     ev2$results)#,
                    # ev3$results)

results_all[results_all$missing_rate == "0.2", "missing_rate"] <- 0.2
results_all[results_all$missing_rate == "0.1", "missing_rate"] <- 0.1
results_all[results_all$missing_rate == "0", "missing_rate"] <- 0
results_all[results_all$sample == "50", "sample"] <- "n = 50"
results_all[results_all$sample == "200", "sample"] <- "n = 200"
results_all[results_all$sample == "500", "sample"] <- "n = 500"

appender <- function(string) TeX(paste("$\\rho_{m} = $", string))

ROC_plot <- ggplot(data = results_all, aes(x = fpr, y = tpr)) +
  geom_line(aes(linetype = method)) +
  scale_linetype_manual(breaks = c("BMGM", "BDgraph", "huge", "MGM"),
                        values=c("solid", "dashed", "dotdash", "dotted")) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color="grey",
              linetype="longdash", size=0.5) +
  labs(x = "False Positive Rate", y = "True Positive Rate", linetype = "") +
  facet_wrap( ~ missing_rate, labeller = as_labeller(appender,
                                                     default = label_parsed))



namefile <- paste0("plot_roc_p_",p,"_n_",n1,"_new.png")
ggsave(namefile, ROC_plot, width = 8, height = 3, device = "png", dpi = 600)


#AUC:

#######################################################

#Corrplot for fixed graph::


p = 14
n1 = n2 = 500
mr1 = 0
mr2 = 0.1

results1 <- read.csv(paste0("C:/Users/User/Documents/Simulation3_fixed_imputation_new_p_",p,"_n_",n1,"_mr_",mr1,".csv"))[,-1]
results2 <- read.csv(paste0("C:/Users/User/Documents/Simulation2_fixed_imputation_new_p_",p,"_n_",n2,"_mr_",mr2,".csv"))[,-1]

#recover the observed graph using the posterior probabilities
#AUC
AUC <- cbind(ev1$auc,
             ev2$auc)#,
             #ev3$auc)

eva1 <- evaluation_table(results1)
eva2 <- evaluation_table(results2)
eva3 <- evaluation_table(results3)

export_tab = t(rbind("AUC" = ev1$auc,eva1[4,],eva1[1:3,]))

xtable(export_tab)

export_tab = t(rbind("AUC" = ev2$auc,eva2[4,],eva2[1:3,]))

xtable(export_tab)

export_tab = t(rbind("AUC" = ev3$auc,eva3[4,],eva3[1:3,]))

xtable(export_tab)



##########################

graphs <- evaluation_fixed(results1)


evaluation_table <- function(file){
  obs_graph = file
  l = p*(p-1)/2

  index <- cbind(seq(1, 5*l, by = l), seq(l, ncol(obs_graph), by = l))

  #HUGE
  obs_graph[,index[1,1]:index[1,2]] <- (abs(obs_graph[,index[1,1]:index[1,2]]) >= 0.5)*1

  #BDgraph

  obs_graph[,index[2,1]:index[2,2]] <- (abs(obs_graph[,index[2,1]:index[2,2]]) >= 0.5)*1

  #MGM

  obs_graph[,index[3,1]:index[3,2]] <- (abs(obs_graph[,index[3,1]:index[3,2]]) > 0)*1

  #BMGM

  #Find the E_FDR

  bmgm_ppd <- abs(obs_graph[,index[4,1]:index[4,2]])
  bfdr = 0.05
  for(i in 1:nrow(bmgm_ppd)){
    post_inclusion <- bmgm_ppd[i,]

    fdr_c <- function(c){
      E_fdr <- sum((1 - post_inclusion)*(post_inclusion > c))/(sum(post_inclusion > c) +
                                                                 rnorm(1, mean = 0, sd = 0.001))
      return(E_fdr)
    }
    pos_c <- seq(0, 1, by = 0.01)
    expected_fdr <- sapply(pos_c, fdr_c)
    pos <- pos_c[expected_fdr < bfdr]
    cutoff <- min(pos)
    obs_graph[i,index[4,1]:index[4,2]] <- (bmgm_ppd[i,] >= cutoff)*1
  }


  #NOW I HAVE THE REPLICATES OF THE GRAPHS!

  #Now let's calculate the measures:

  mcc = spec = sens = f1 = tpr = fpr = matrix(nrow = nrow(obs_graph), ncol = 4)

  for(i in 1:nrow(obs_graph)){
    adj <- diag(p) - diag(p)

    adj_huge = adj
    adj_huge[upper.tri(adj_huge)] <- unlist(obs_graph[i, index[1,1]:index[1,2]])
    adj_huge <- adj_huge + t(adj_huge)

    adj_bd = adj
    adj_bd[upper.tri(adj_bd)] <- unlist(obs_graph[i, index[2,1]:index[2,2]])
    adj_bd <- adj_bd + t(adj_bd)

    adj_mgm = adj
    adj_mgm[upper.tri(adj_mgm)] <- unlist(obs_graph[i, index[3,1]:index[3,2]])
    adj_mgm <- adj_mgm + t(adj_mgm)

    adj_bmgm = adj
    adj_bmgm[upper.tri(adj_bmgm)] <- unlist(obs_graph[i, index[4,1]:index[4,2]])
    adj_bmgm <- adj_bmgm + t(adj_bmgm)

    true <- adj
    true[upper.tri(true)] <-  unlist(obs_graph[i, index[5,1]:index[5,2]])
    true <- true + t(true)

    comparison <- BDgraph::compare(list(adj_huge, adj_bd, adj_mgm, adj_bmgm), true)[,-1]

    mcc[i,] <- comparison["MCC",]
    spec[i,] <- comparison["Specificity",]
    sens[i,] <- comparison["Sensitivity",]
    f1[i,] <- comparison["F1-score",]
    tpr[i,] <- comparison["True Positive", ] / (comparison["True Positive", ] + comparison["False Negative", ])
    fpr[i,] <- comparison["False Positive", ] / (comparison["False Positive", ] + comparison["True Negative", ])
  }

  summ <- rbind(apply(mcc, 2, mean),
                apply(spec, 2, mean),
                apply(sens, 2, mean),
                apply(f1, 2, mean),
                apply(tpr, 2, mean),
                apply(fpr, 2, mean))

  summ_sd <- rbind(apply(mcc, 2, sd),
                apply(spec, 2, sd),
                apply(sens, 2, sd),
                apply(f1, 2, sd),
                apply(tpr, 2, sd),
                apply(fpr, 2, sd))

  colnames(summ) = colnames(summ_sd) <- c("huge", "BDgraph", "MGM", "BMBM")

  rownames(summ) = rownames(summ_sd) <- c("MCC", "Specificity", "Sensitivity", "F1", "TPR", "FPR")

  summ[] <- paste0(round(summ, 3), " (", round(summ_sd, 3), ")")

  return(summ)
}

evaluation <- function(file, n, p, mr){
  l = p*(p-1)/2
  results <- file

  mcc_final <- c()

  results.huge <- results[,1:l]
  true <- results[, (ncol(results) - l + 1):ncol(results)]

  #for each iteration, save the threshold, fdr, and tpr
  rep <- nrow(results)
  th <- seq(-0.1, 1, by = 0.01)

  tpr <- matrix(nrow = rep, ncol = length(th))
  fpr <- matrix(nrow = rep, ncol = length(th))
  mcc <- matrix(nrow = rep, ncol = length(th))

  for(i in 1:rep){
    true_i <- diag(p) - diag(p)
    true_i[upper.tri(true_i)] <- as.numeric(true[i,])
    true_i <- true_i + t(true_i)
    true_i <- abs(true_i)

    for(t in 1:length(th)){
      adj <- diag(p) - diag(p)
      adj[upper.tri(adj)] <- (abs(results.huge[i,]) > th[t])*1
      adj <- adj + t(adj)
      comp <- BDgraph::compare(adj, true_i)
      tpr[i,t] <- comp["True Positive", 2] / (comp["True Positive", 2] + comp["False Negative", 2])
      fpr[i,t] <- comp["False Positive", 2] / (comp["False Positive", 2] + comp["True Negative", 2])
      mcc[i,t] <- comp["MCC", 2]
    }
  }

  mcc_final[1] <- max(colMeans(mcc))
  plot.tpr.huge <- colMeans(tpr)
  plot.fpr.huge <- colMeans(fpr)

  auc <- rep(0, nrow(fpr))
  for(i in 1:length(auc)){
    auc[i] <- flux::auc(fpr[i,], tpr[i,])
  }
  auc_huge <- c(mean(auc), sd(auc))

  ########################### BDGRAPH ################

  results.bd <- results[,(l+1):(2*l)]
  true <- results[, (ncol(results) - l + 1):ncol(results)]

  #for each iteration, save the threshold, fdr, and tpr
  rep <- nrow(results)
  th <- sort(seq(-0.1, 1, by = 0.01), decreasing = T)

  tpr <- matrix(nrow = rep, ncol = length(th))
  fpr <- matrix(nrow = rep, ncol = length(th))
  mcc <- matrix(nrow = rep, ncol = length(th))


  for(i in 1:rep){
    true_i <- diag(p) - diag(p)
    true_i[upper.tri(true_i)] <- as.numeric(true[i,])
    true_i <- true_i + t(true_i)
    true_i <- abs(true_i)

    for(t in 1:length(th)){
      adj <- diag(p) - diag(p)
      adj[upper.tri(adj)] <- (abs(results.bd[i,]) > th[t])*1
      adj <- adj + t(adj)
      comp <- BDgraph::compare(adj, true_i)
      tpr[i,t] <- comp["True Positive", 2] / (comp["True Positive", 2] + comp["False Negative", 2])
      fpr[i,t] <- comp["False Positive", 2] / (comp["False Positive", 2] + comp["True Negative", 2])
      mcc[i,t] <- comp["MCC", 2]
    }
  }

  mcc_final[2] <- max(colMeans(mcc))
  plot.tpr.bd <- colMeans(tpr)
  plot.fpr.bd <- colMeans(fpr)
  auc <- rep(0, nrow(fpr))
  for(i in 1:length(auc)){
    auc[i] <- flux::auc(fpr[i,], tpr[i,])
  }
  auc_bd <- c(mean(auc), sd(auc))

  ########################### MGM ################

  results.mgm <- results[,(2*l+1):(3*l)]
  true <- results[, (ncol(results) - l + 1):ncol(results)]


  #for each iteration, save the threshold, fdr, and tpr
  rep <- nrow(results)
  th <- seq(-0.1, 1, by = 0.001)

  tpr <- matrix(nrow = rep, ncol = length(th))
  fpr <- matrix(nrow = rep, ncol = length(th))
  mcc <- matrix(nrow = rep, ncol = length(th))


  for(i in 1:rep){
    true_i <- diag(p) - diag(p)
    true_i[upper.tri(true_i)] <- as.numeric(true[i,])
    true_i <- true_i + t(true_i)
    true_i <- abs(true_i)

    for(t in 1:length(th)){
      adj <- diag(p) - diag(p)
      adj[upper.tri(adj)] <- (abs(results.mgm[i,]) > th[t])*1
      adj <- adj + t(adj)
      comp <- BDgraph::compare(adj, true_i)
      tpr[i,t] <- comp["True Positive", 2] / (comp["True Positive", 2] + comp["False Negative", 2])
      fpr[i,t] <- comp["False Positive", 2] / (comp["False Positive", 2] + comp["True Negative", 2])
      mcc[i,t] <- comp["MCC", 2]
    }
  }

  mcc_final[3] <- max(colMeans(mcc))
  plot.tpr.mgm <- colMeans(tpr)
  plot.fpr.mgm <- colMeans(fpr)
  auc <- rep(0, nrow(fpr))
  for(i in 1:length(auc)){
    auc[i] <- flux::auc(fpr[i,], tpr[i,])
  }
  auc_mgm <- c(mean(auc), sd(auc))

  ######################### OUR #######################

  results.our <- results[,(3*l+1):(4*l)]
  true <- results[, (ncol(results) - l + 1):ncol(results)]

  #for each iteration, save the threshold, fdr, and tpr
  rep <- nrow(results)
  th <- seq(-0.1, 1, by = 0.01)

  tpr <- matrix(nrow = rep, ncol = length(th))
  fpr <- matrix(nrow = rep, ncol = length(th))
  mcc <- matrix(nrow = rep, ncol = length(th))

  for(i in 1:rep){
    true_i <- diag(p) - diag(p)
    true_i[upper.tri(true_i)] <- as.numeric(true[i,])
    true_i <- true_i + t(true_i)
    true_i <- abs(true_i)

    for(t in 1:length(th)){
      adj <- diag(p) - diag(p)
      adj[upper.tri(adj)] <- c((abs(results.our[i,]) > th[t])*1)
      adj <- adj + t(adj)
      comp <- BDgraph::compare(adj, true_i)
      tpr[i,t] <- comp["True Positive", 2] / (comp["True Positive", 2] + comp["False Negative", 2])
      fpr[i,t] <- comp["False Positive", 2] / (comp["False Positive", 2] + comp["True Negative", 2])
      mcc[i,t] <- comp["MCC", 2]
    }
  }

  #ROC curve:

  plot.tpr.our <- colMeans(tpr)
  plot.fpr.our <- colMeans(fpr)
  mcc_final[4] <- max(colMeans(mcc))
  auc <- rep(0, nrow(fpr))
  for(i in 1:length(auc)){
    auc[i] <- flux::auc(fpr[i,], tpr[i,])
  }
  auc_bmgm <- c(mean(auc), sd(auc))


  all_tpr <- c(plot.tpr.huge, plot.tpr.bd, plot.tpr.mgm, plot.tpr.our)
  all_fpr <- c(plot.fpr.huge, plot.fpr.bd, plot.fpr.mgm, plot.fpr.our)


  auc_all <- round(rbind(auc_huge, auc_bd, auc_mgm, auc_bmgm), 3)
  auc_all <- t(as.matrix(paste0(auc_all[,1], " (", auc_all[,2], ")")))
  colnames(auc_all) <- c("huge", "BDgraph", "MGM", "BMGM")

  method <- c(rep("huge", length(plot.fpr.huge)), rep("BDgraph", length(plot.fpr.bd)),
              rep("MGM", length(plot.fpr.mgm)), rep("BMGM", length(plot.fpr.our)))


  sample_size <- rep(n, length(method))
  missing_rate <- rep(mr, length(method))

  #results_all <- data.frame(tpr = all_tpr, fpr = all_fpr, method = method, sample = sample_size)
  results_all <- data.frame(tpr = all_tpr, fpr = all_fpr, method = method,
                            sample = sample_size, missing_rate = missing_rate)

  return(list(results = results_all, auc = auc_all))
}

evaluation_fixed <- function(file){
  obs_graph = file[,c(1:91,547:728, 820:910, 365:455)]
  l = p*(p-1)/2

  index <- cbind(seq(1, 5*l, by = l), seq(l, ncol(obs_graph), by = l))

  #HUGE
  obs_graph[,index[1,1]:index[1,2]] <- (abs(obs_graph[,index[1,1]:index[1,2]]) >= 0.5)*1

  #BDgraph

  #obs_graph[,index[2,1]:index[2,2]] <- (abs(obs_graph[,index[2,1]:index[2,2]]) >= 0.5)*1

  #MGM

  #obs_graph[,index[3,1]:index[3,2]] <- (abs(obs_graph[,index[3,1]:index[3,2]]) > 0)*1

  #BMGM

  #Find the E_FDR

  bmgm_ppd <- abs(obs_graph[,index[4,1]:index[4,2]])
  bfdr = 0.01
  for(i in 1:nrow(bmgm_ppd)){
    post_inclusion <- bmgm_ppd[i,]

    fdr_c <- function(c){
      E_fdr <- sum((1 - post_inclusion)*(post_inclusion > c))/(sum(post_inclusion > c) +
                                                                 rnorm(1, mean = 0, sd = 0.001))
      return(E_fdr)
    }
    pos_c <- seq(0, 1, by = 0.01)
    expected_fdr <- sapply(pos_c, fdr_c)
    pos <- pos_c[expected_fdr < bfdr]
    cutoff <- min(pos)
    obs_graph[i,index[4,1]:index[4,2]] <- (bmgm_ppd[i,] > cutoff)*1
  }

  return(obs_graph)
}

graphs <- evaluation_fixed(results1)

graphs_mean <- colMeans(graphs)

l = p*(p-1)/2

index <- cbind(seq(1, 5*l, by = l), seq(l, length(graphs_mean), by = l))

adj <- diag(p) - diag(p)

adj_huge = adj
adj_huge[upper.tri(adj_huge)] <- unlist(graphs_mean[index[1,1]:index[1,2]])
adj_huge <- adj_huge + t(adj_huge)

adj_bd = adj
adj_bd[upper.tri(adj_bd)] <- unlist(graphs_mean[index[2,1]:index[2,2]])
adj_bd <- adj_bd + t(adj_bd)

adj_mgm = adj
adj_mgm[upper.tri(adj_mgm)] <- unlist(graphs_mean[index[3,1]:index[3,2]])
adj_mgm <- adj_mgm + t(adj_mgm)

adj_bmgm = adj
adj_bmgm[upper.tri(adj_bmgm)] <-  unlist(graphs_mean[index[4,1]:index[4,2]])
adj_bmgm <- adj_bmgm + t(adj_bmgm)

true <- adj
true[upper.tri(true)] <-  unlist(graphs_mean[index[5,1]:index[5,2]])
true <- true + t(true)

#######

corrplot(true, method = "square")

g1 <- ggcorrplot(true,
                 ggtheme = ggplot2::theme_void,
                 show.legend = T) +
  scale_fill_gradient2(breaks=c(0, 0.25, 0.5, 0.75, 1), limit=c(0, 1), low = "white", high = "black") +
  labs(fill = "")

leg <- get_legend(g1)

g1 <- ggcorrplot(adj_huge, title = "Huge",
                 ggtheme = ggplot2::theme_void,
                 show.legend = F) +
  scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1), low = "white", high = "black") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

g2 <- ggcorrplot(adj_bd, title = "BDgraph",
                 ggtheme = ggplot2::theme_void,
                 show.legend = F) +
  scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1), low = "white", high = "black") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

g3 <- ggcorrplot(adj_mgm, title = "MGM",
                 ggtheme = ggplot2::theme_void,
                 show.legend = F) +
  scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1), low = "white", high = "black") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())


g4 <- ggcorrplot(adj_bmgm, title = "BMGM",
                 ggtheme = ggplot2::theme_void,
                 show.legend = F) +
  scale_fill_gradient2(breaks=c(0, 0.25, 0.5, 0.75, 1), limit=c(0, 1), low = "white", high = "black") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

g0 <- ggcorrplot(true, title = "Observed Graph",
                 ggtheme = ggplot2::theme_void,
                 show.legend = F) +
  labs(x = "Node", y = "Node") +
  scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1), low = "white", high = "black") +
  theme(plot.title = element_text(hjust = 0.5))


PPI1 <- plot_grid(NULL,g0, NULL, g4, g2, g1, g3, leg, nrow = 1,
                  rel_widths = c(0.2,1.1,0.2,1,1,1,1,0.5))
PPI1

namefile <- paste0("PPI_p14_mr0_n500.png")
ggsave(namefile, PPI, width = 10, height = 3, device = "png", dpi = 600)

namefile <- paste0("PPI_p14_mr0.1_n500.png")
ggsave(namefile, PPI2, width = 10, height = 3, device = "png", dpi = 600)

###############################################################################


# FOR ONE REPLICATEp = 14


library(ggcorrplot)
library(cowplot)
library(corrplot)


PPI_plot <- function(file, repl = 2){
  results1 = file
  graphs = evaluation_fixed(results1)
  # Number of edges
  p <- 14
  l <- choose(p, 2)

  # Recover the hard selections (each block of l columns)
  hard_huge   <- as.matrix(graphs[, (0*l + 1):(1*l)])
  hard_bd   <- as.matrix(graphs[, (1*l + 1):(2*l)])
  hard_mgm     <- as.matrix(graphs[, (2*l + 1):(3*l)])
  hard_our    <- as.matrix(graphs[, (3*l + 1):(4*l)])
  # skip (8–9)*l since it's repeated BDgraph
  true_graph    <- as.matrix(graphs[, (4*l + 1):(5*l)] > 0)*1

  library(ggcorrplot)
  library(cowplot)
  library(corrplot)

  # Replicate index
  adj <- matrix(0, nrow = p, ncol = p)
  type <- c("B", "C", "Z", "C", "Z", "C", "Z", "B", "C", "Z", "C", "Z", "C", "Z")


  make_adj <- function(vec, truev) {
    # Mark false negatives with -1
    vec[vec == 1 & truev == 0] <- -1

    # Create adjacency matrix
    mat <- matrix(0, p, p)
    mat[upper.tri(mat)] <- vec
    mat <- mat + t(mat)
    diag(mat) <- 0
    return(mat)
  }

  adj_huge  <- make_adj(hard_huge[repl, ], true_graph[repl,])
  adj_bd    <- make_adj(hard_bd[repl, ], true_graph[repl,])
  adj_mgm   <- make_adj(hard_mgm[repl, ], true_graph[repl,])
  adj_bmgm  <- make_adj(hard_our[repl, ], true_graph[repl,])
  true      <- make_adj(true_graph[repl, ], true_graph[repl,])


  g1 <- ggcorrplot(true,
                   ggtheme = ggplot2::theme_void,
                   show.legend = T) +
    scale_fill_gradient2(breaks=c(0, 0.25, 0.5, 0.75, 1), limit=c(-1, 1), low = "firebrick", high = "black") +
    labs(fill = "")

  leg <- get_legend(g1)

  g1 <- ggcorrplot(adj_huge, title = "Huge",
                   ggtheme = ggplot2::theme_void,
                   show.legend = F) +
    scale_fill_gradient2(breaks=c(0, 1), limit=c(-1, 1), low = "firebrick", high = "black") +
    scale_x_continuous(breaks = 1:p, labels = type) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0, vjust = 0.5, size = 8, hjust = 0.5),
          axis.text.y=element_blank())

  g2 <- ggcorrplot(adj_bd, title = "BDgraph",
                   ggtheme = ggplot2::theme_void,
                   show.legend = F) +
    scale_fill_gradient2(breaks=c(0, 1), limit=c(-1, 1), low = "firebrick", high = "black") +
    scale_x_continuous(breaks = 1:p, labels = type) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0, vjust = 0.5, size = 8, hjust = 0.5),
          axis.text.y=element_blank())

  g3 <- ggcorrplot(adj_mgm, title = "MGM",
                   ggtheme = ggplot2::theme_void,
                   show.legend = F) +
    scale_fill_gradient2(breaks=c(0, 1), limit=c(-1, 1), low = "firebrick", high = "black") +
    scale_x_continuous(breaks = 1:p, labels = type) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0, vjust = 0, size = 8, hjust = 0.5),
          axis.text.y=element_blank())


  g4 <- ggcorrplot(adj_bmgm, title = "BMGM",
                   ggtheme = ggplot2::theme_void,
                   show.legend = F) +
    scale_fill_gradient2(breaks=c(0, 1), limit=c(-1, 1), low = "firebrick", high = "black") +
    scale_x_continuous(breaks = 1:p, labels = type) +
    scale_y_continuous(breaks = 1:p, labels = type) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0, vjust = 0.5, size = 8, hjust = 0.5),
          axis.text.y= element_text(size = 8))

  g0 <- ggcorrplot(true, title = "Observed Graph",
                   ggtheme = ggplot2::theme_void,
                   show.legend = F) +
    labs(x = "Node Type", y = "Node Type") +
    scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1), low = "white", high = "black") +
    scale_x_continuous(breaks = 1:p, labels = type) +
    scale_y_continuous(breaks = 1:p, labels = type) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0, vjust = 0.5, size = 8, hjust = 0.5),
          axis.text.y = element_text(size = 8))



  PPI1 <- plot_grid(
    NULL, g0, NULL, g4, g2, g1, g3,
    nrow = 1,
    rel_widths = c(0.2, 1.1, 0.2, 1, 1, 1, 1, 0.5)
  )
  return(PPI1)
}

l = p*(p-1)/2
mr1 = 0
mr2 = 0.1

n1 = 200
n2 = 200

#results1 <- read.csv(paste0("C:/Users/User/Documents/Simulation_imputation_p_",p,"_n_",n1,"_mr_",mr1,".csv"))[,-1]
results1 <- read.csv(paste0("C:/Users/User/Documents/Simulation2_fixed_imputation_new_p_",p,"_n_",n1,"_mr_",mr1,".csv"))[,-1]
results2 <- read.csv(paste0("C:/Users/User/Documents/Simulation2_fixed_imputation_new_p_",p,"_n_",n2,"_mr_",mr2,".csv"))[,-1]

n1 = 500
n2 = 500

results3 <- read.csv(paste0("C:/Users/User/Documents/Simulation2_fixed_imputation_new_p_",p,"_n_",n1,"_mr_",mr1,".csv"))[,-1]
results4 <- read.csv(paste0("C:/Users/User/Documents/Simulation2_fixed_imputation_new_p_",p,"_n_",n2,"_mr_",mr2,".csv"))[,-1]

PPIs <- list()
PPIs[[1]] <- PPI_plot(results1)
PPIs[[2]] <- PPI_plot(results2)
PPIs[[3]] <- PPI_plot(results3)
PPIs[[4]] <- PPI_plot(results4)

library(ggplot2)
library(cowplot)

# Create dummy data for the legend
# Data for each legend entry
legend_tp <- data.frame(x = 1, y = 1, label = "True Positive", color = "black")
legend_fn <- data.frame(x = 1, y = 1, label = "False Positive", color = "firebrick")

# Each legend item as a separate ggplot
legend_tp_plot <- ggplot(legend_tp, aes(x, y, fill = color)) +
  geom_tile(width = 0.1, height = 0.4) +
  geom_text(aes(label = label), hjust = -0.2, size = 4) +
  scale_fill_identity() +
  theme_void() + xlim(0.5, 2.5) + ylim(0.5, 1.5)

legend_fn_plot <- ggplot(legend_fn, aes(x, y, fill = color)) +
  geom_tile(width = 0.1, height = 0.4) +
  geom_text(aes(label = label), hjust = -0.2, size = 4) +
  scale_fill_identity() +
  theme_void() + xlim(0.5, 2.5) + ylim(0.5, 1.5)

legend_combined <- plot_grid(NULL,NULL, legend_tp_plot, legend_fn_plot, NULL, NULL,
                             nrow = 1,
                             rel_widths = c(1, 1, 1.5, 1.5, 1, 1))

final_plot <- plot_grid(PPIs[[1]], PPIs[[2]], PPIs[[3]], PPIs[[4]], legend_combined, ncol = 1, rel_heights = c(1,1,1,1,0.2))
final_plot

namefile <- paste0("PPI_onereplicate_p14.png")
ggsave(namefile, final_plot, width = 10, height = 8, device = "png", dpi = 600)



###################################
#METRICS

calculate_metrics_all <- function(file){
  results1 = file
  graphs = evaluation_fixed(results1)

  p <- 14
  l <- choose(p, 2)

  hard_huge <- as.matrix(graphs[, (0*l + 1):(1*l)])
  hard_bd <- as.matrix(graphs[, (1*l + 1):(2*l)])
  hard_mgm <- as.matrix(graphs[, (2*l + 1):(3*l)])
  hard_our <- as.matrix(graphs[, (3*l + 1):(4*l)])
  true_graph <- as.matrix(graphs[, (4*l + 1):(5*l)] > 0) * 1

  make_adj <- function(vec) {
    mat <- matrix(0, p, p)
    mat[upper.tri(mat)] <- vec
    mat <- mat + t(mat)
    diag(mat) <- 0
    return(mat)
  }

  library(BDgraph)

  metrics <- data.frame()

  for (repl in 1:nrow(hard_huge)) {
    adj_huge  <- make_adj(hard_huge[repl, ])
    adj_bd    <- make_adj(hard_bd[repl, ])
    adj_mgm   <- make_adj(hard_mgm[repl, ])
    adj_bmgm  <- make_adj(hard_our[repl, ])
    adj_true  <- make_adj(true_graph[repl, ])

    cmp <- list(
      Huge = BDgraph::compare(adj_true, adj_huge),
      BDgraph = BDgraph::compare(adj_true, adj_bd),
      MGM = BDgraph::compare(adj_true, adj_mgm),
      BMGM = BDgraph::compare(adj_true, adj_bmgm)
    )

    for (method in names(cmp)) {
      res <- cmp[[method]]
      metrics <- rbind(metrics, c(Method = method, res[,2]))
    }
  }
  colnames(metrics) <- c("Method", "TP", "TN", "FP", "FN", "F1-score", "Specificity", "Sensitivity", "MCC")
  metrics[,-1] <- sapply(metrics[,-1], as.numeric)
  aggregate(. ~ Method, data = metrics, FUN = function(x) sprintf("%.3f (%.3f)", mean(x), sd(x)))
}

calculate_metrics_all(results1)

calculate_metrics_all(results2)


################################################################################

# sorting the plot:

library(ggcorrplot)
library(cowplot)
library(corrplot)

# Setup
p <- 14
l <- choose(p, 2)

PPI_plot <- function(file, repl = 2){
# Load and extract matrices
  results1 = file
  graphs = evaluation_fixed(results1)

  hard_huge   <- as.matrix(graphs[, (0*l + 1):(1*l)])
  hard_bd     <- as.matrix(graphs[, (1*l + 1):(2*l)])
  hard_mgm    <- as.matrix(graphs[, (2*l + 1):(3*l)])
  hard_our    <- as.matrix(graphs[, (3*l + 1):(4*l)])
  true_graph  <- as.matrix(graphs[, (4*l + 1):(5*l)] > 0)*1

  # Original type vector
  type <- c("B", "C", "Z", "C", "Z", "C", "Z", "B", "C", "Z", "C", "Z", "C", "Z")

  # Reorder indices: Binary → Continuous → Zero-inflated
  new_order <- order(factor(type, levels = c("B", "C", "Z")))
  type_ordered <- type[new_order]

  # Function to convert edge vector to adjacency matrix with TP/FP marking
  make_adj <- function(vec, truev) {
    vec[vec == 1 & truev == 0] <- -1
    mat <- matrix(0, p, p)
    mat[upper.tri(mat)] <- vec
    mat <- mat + t(mat)
    diag(mat) <- 0
    return(mat)
  }

  # Build and reorder adjacency matrices
  adj_huge  <- make_adj(hard_huge[repl, ], true_graph[repl,])
  adj_bd    <- make_adj(hard_bd[repl, ], true_graph[repl,])
  adj_mgm   <- make_adj(hard_mgm[repl, ], true_graph[repl,])
  adj_bmgm  <- make_adj(hard_our[repl, ], true_graph[repl,])
  true      <- make_adj(true_graph[repl, ], true_graph[repl,])

  # Reorder adjacency matrices and true matrix
  reorder_adj <- function(mat, order) mat[order, order]
  adj_huge  <- reorder_adj(adj_huge, new_order)
  adj_bd    <- reorder_adj(adj_bd, new_order)
  adj_mgm   <- reorder_adj(adj_mgm, new_order)
  adj_bmgm  <- reorder_adj(adj_bmgm, new_order)
  true      <- reorder_adj(true, new_order)

  # Common plotting theme
  plot_theme <- theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 0, vjust = 0.5, size = 8, hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )

  # Plot for each method
  plot_corr <- function(mat, title, show_y = FALSE) {
    ggcorrplot(mat, title = title, ggtheme = ggplot2::theme_void, show.legend = F) +
      scale_fill_gradient2(breaks=c(-1, 0, 1), limit=c(-1, 1), low = "firebrick", mid = "white", high = "black") +
      scale_x_continuous(breaks = 1:p, labels = type_ordered) +
      scale_y_continuous(breaks = if (show_y) 1:p else NULL, labels = if (show_y) type_ordered else NULL) +
      plot_theme +
      guides(fill = "none")
  }

  # Create individual plots
  g0 <- ggcorrplot(true, title = "Observed Graph",
                   ggtheme = theme_void(), show.legend = F) +
    labs(x = "Node Type", y = "Node Type") +
    scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1), low = "white", high = "black") +
    scale_x_continuous(breaks = 1:p, labels = type_ordered) +
    scale_y_continuous(breaks = 1:p, labels = type_ordered) +
    plot_theme +
    guides(fill = "none")

  g1 <- plot_corr(adj_huge, "Huge")
  g2 <- plot_corr(adj_bd, "BDgraph")
  g3 <- plot_corr(adj_mgm, "MGM")
  g4 <- plot_corr(adj_bmgm, "BMGM", show_y = TRUE)

  # Combine plots
  PPI1 <- plot_grid(
    NULL, g0, NULL, g4, g2, g1, g3,
    nrow = 1,
    rel_widths = c(0.2, 1.1, 0.2, 1, 1, 1, 1, 0.5)
  )

  # Return or print the final plot grid
  PPI1
}


l = p*(p-1)/2
mr1 = 0
mr2 = 0.1

n1 = 200
n2 = 200

#results1 <- read.csv(paste0("C:/Users/User/Documents/Simulation_imputation_p_",p,"_n_",n1,"_mr_",mr1,".csv"))[,-1]
results1 <- read.csv(paste0("C:/Users/User/Documents/Simulation2_fixed_imputation_new_p_",p,"_n_",n1,"_mr_",mr1,".csv"))[,-1]
results2 <- read.csv(paste0("C:/Users/User/Documents/Simulation2_fixed_imputation_new_p_",p,"_n_",n2,"_mr_",mr2,".csv"))[,-1]

n1 = 500
n2 = 500

results3 <- read.csv(paste0("C:/Users/User/Documents/Simulation2_fixed_imputation_new_p_",p,"_n_",n1,"_mr_",mr1,".csv"))[,-1]
results4 <- read.csv(paste0("C:/Users/User/Documents/Simulation2_fixed_imputation_new_p_",p,"_n_",n2,"_mr_",mr2,".csv"))[,-1]

PPIs <- list()
PPIs[[1]] <- PPI_plot(results1)
PPIs[[2]] <- PPI_plot(results2)
PPIs[[3]] <- PPI_plot(results3)
PPIs[[4]] <- PPI_plot(results4)

legend_tp_plot <- ggplot(legend_tp, aes(x, y, fill = color)) +
  geom_tile(width = 0.1, height = 0.4) +
  geom_text(aes(label = label), hjust = -0.2, size = 4) +
  scale_fill_identity() +
  theme_void() + xlim(0.5, 2.5) + ylim(0.5, 1.5)

legend_fn_plot <- ggplot(legend_fn, aes(x, y, fill = color)) +
  geom_tile(width = 0.1, height = 0.4) +
  geom_text(aes(label = label), hjust = -0.2, size = 4) +
  scale_fill_identity() +
  theme_void() + xlim(0.5, 2.5) + ylim(0.5, 1.5)

legend_combined <- plot_grid(NULL,NULL, legend_tp_plot, legend_fn_plot, NULL, NULL,
                             nrow = 1,
                             rel_widths = c(1, 1, 1.5, 1.5, 1, 1))

final_plot <- plot_grid(PPIs[[1]], PPIs[[2]], PPIs[[3]], PPIs[[4]], legend_combined, ncol = 1, rel_heights = c(1,1,1,1,0.2))
final_plot

namefile <- paste0("PPI_onereplicate_arranged_p14.png")
ggsave(namefile, final_plot, width = 10, height = 8, device = "png", dpi = 600)


