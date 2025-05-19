library(dplyr)
library(ggplot2)

# Load
#ppi_df <- read.csv("C:/Users/User/Documents/GitHub/BMGM/simulations/Simulation_bdgraph_p_10_N_200_mr_0.csv")
#ppi_df <- read.csv("C:/Users/User/Documents/GitHub/BMGM/simulations/Simulation_fixed_imputation_mgm2_p_10_N_200_mr_0.1.csv")

ppi_df <- read.csv("C:/Users/User/Documents/GitHub/BMGM/simulations/Simulation_null_p_10_N_500_mr_0.2.csv")

ppi_df <- ppi_df[, -1]

# Setup
ncols <- choose(10,2)  # because p = 10

# Now assign properly
soft_huge   <- as.matrix(ppi_df[, 1:ncols])
soft_bd     <- as.matrix(ppi_df[, (ncols+1):(2*ncols)])
soft_mgm    <- as.matrix(ppi_df[, (2*ncols+1):(3*ncols)])
soft_our1   <- as.matrix(ppi_df[, (3*ncols+1):(4*ncols)])
soft_our2   <- as.matrix(ppi_df[, (4*ncols+1):(5*ncols)])
true_graph  <- as.matrix(ppi_df[, (5*ncols+1):(6*ncols)])

# Hard selections
hard_huge   <- as.matrix(ppi_df[, (6*ncols+1):(7*ncols)])
hard_bd     <- as.matrix(ppi_df[, (7*ncols+1):(8*ncols)])
hard_mgm    <- as.matrix(ppi_df[, (8*ncols+1):(9*ncols)])
# (9–10 block is the repeated bdgraph, ignore)
hard_our1   <- as.matrix(ppi_df[, (10*ncols+1):(11*ncols)])
hard_our2   <- as.matrix(ppi_df[, (11*ncols+1):(12*ncols)])

compute_roc_replicate <- function(pred_matrix, truth_matrix, thresholds = seq(0, 1, 0.01)) {
  nrep <- nrow(pred_matrix)
  all_roc <- vector("list", nrep)

  for (i in 1:nrep) {
    pred <- pred_matrix[i, ]
    truth <- truth_matrix[i, ]

    roc_i <- lapply(thresholds, function(t) {
      pred_bin <- ifelse(pred >= t, 1, 0)
      tp <- sum(pred_bin == 1 & truth == 1)
      fp <- sum(pred_bin == 1 & truth == 0)
      fn <- sum(pred_bin == 0 & truth == 1)
      tn <- sum(pred_bin == 0 & truth == 0)
      tpr <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
      fpr <- ifelse(fp + tn == 0, 0, fp / (fp + tn))
      data.frame(fpr = fpr, tpr = tpr)
    }) %>% bind_rows()

    all_roc[[i]] <- roc_i
  }

  return(all_roc)
}

average_roc <- function(all_roc_list) {
  thresholds_len <- nrow(all_roc_list[[1]])
  nrep <- length(all_roc_list)

  mean_tpr <- c(sapply(1:thresholds_len, function(j) mean(sapply(all_roc_list, function(x) x$tpr[j]))),0)
  mean_fpr <- c(sapply(1:thresholds_len, function(j) mean(sapply(all_roc_list, function(x) x$fpr[j]))),0)

  data.frame(fpr = mean_fpr, tpr = mean_tpr)
}

# Compute replicate-by-replicatehttp://127.0.0.1:40889/graphics/plot_zoom_png?width=2698&height=1093
roc_huge_reps  <- compute_roc_replicate(abs(soft_huge),  true_graph)
roc_bd_reps    <- compute_roc_replicate(abs(soft_bd),    true_graph)
roc_mgm_reps   <- compute_roc_replicate(abs(soft_mgm),   true_graph)
roc_our1_reps  <- compute_roc_replicate(abs(soft_our1),  true_graph)
roc_our2_reps  <- compute_roc_replicate(abs(soft_our2),  true_graph)

# Average across replicates
roc_huge  <- average_roc(roc_huge_reps) |> mutate(method = "Huge")
roc_bd    <- average_roc(roc_bd_reps)   |> mutate(method = "BDgraph")
roc_mgm   <- average_roc(roc_mgm_reps)  |> mutate(method = "MGM")
roc_our1  <- average_roc(roc_our1_reps) |> mutate(method = "BMGM_v1")
roc_our2  <- average_roc(roc_our2_reps) |> mutate(method = "BMGM_v2")

roc_huge <- roc_huge[order(roc_huge$fpr), ]
roc_bd    <- roc_bd[order(roc_bd$fpr), ]
roc_mgm   <- roc_mgm[order(roc_mgm$fpr), ]
roc_our1  <- roc_our1[order(roc_our1$fpr), ]
roc_our2  <- roc_our2[order(roc_our2$fpr), ]


roc_data <- bind_rows(roc_huge, roc_bd, roc_mgm, roc_our1, roc_our2)

ROC_plot <- ggplot(data = roc_data, aes(x = fpr, y = tpr)) +
  geom_line(aes(linetype = method)) +
  scale_linetype_manual(breaks = c("BMGM_v1", "BMGM_v2", "BDgraph", "Huge", "MGM"),
                        values=c("solid","dotted", "dashed", "dotdash", "dotted")) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color="grey",
              linetype="longdash", size=0.5) +
  labs(x = "False Positive Rate", y = "True Positive Rate", linetype = "")

ROC_plot

compute_auc <- function(roc_df) {
  roc_df <- roc_df[order(roc_df$fpr), ]  # Ensure sorted by FPR
  x <- roc_df$fpr
  y <- roc_df$tpr
  auc <- sum(diff(x) * head(y, -1))  # Trapezoidal integration
  return(auc)
}

# AUCs for each method
auc_huge  <- compute_auc(roc_huge)
auc_bd    <- compute_auc(roc_bd)
auc_mgm   <- compute_auc(roc_mgm)
auc_our1  <- compute_auc(roc_our1)
auc_our2  <- compute_auc(roc_our2)

# Summary
AUC <- round(c(Huge = auc_huge, BDgraph = auc_bd, MGM = auc_mgm, BMGM_v1 = auc_our1, BMGM_v2 = auc_our2), 3)
print(AUC)

# Using the hard selections
hard_huge   <- as.matrix(ppi_df[, (6*ncols+1):(7*ncols)])
hard_bd     <- as.matrix(ppi_df[, (7*ncols+1):(8*ncols)])
hard_mgm    <- as.matrix(ppi_df[, (8*ncols+1):(9*ncols)])
# careful: (9*ncols+1):(10*ncols) is repeated bdgraph -- ignore it
hard_our1   <- as.matrix(ppi_df[, (10*ncols+1):(11*ncols)])
hard_our2   <- as.matrix(ppi_df[, (11*ncols+1):(12*ncols)])

compute_metrics_extended <- function(pred_matrix, truth_matrix) {
  nrep <- nrow(pred_matrix)

  results <- matrix(0, nrow = nrep, ncol = 6)  # TPR, FPR, Precision, F1, MCC

  for (i in 1:nrep) {
    pred <- pred_matrix[i, ]
    truth <- truth_matrix[i, ]

    tp <- sum(pred == 1 & truth == 1)
    fp <- sum(pred == 1 & truth == 0)
    fn <- sum(pred == 0 & truth == 1)
    tn <- sum(pred == 0 & truth == 0)

    tpr <- ifelse(tp + fn == 0, 0, tp / (tp + fn))  # Recall
    fpr <- ifelse(fp + tn == 0, 0, fp / (fp + tn))
    prec <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    f1 <- ifelse(prec + tpr == 0, 0, 2 * prec * tpr / (prec + tpr))

    mcc_denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc <- ifelse(mcc_denominator == 0, 0, (tp * tn - fp * fn) / mcc_denominator)

    results[i, ] <- c(tpr, fpr, prec, f1, mcc, tn)
  }

  colnames(results) <- c("TPR", "FPR", "Precision", "F1", "MCC", "TN")
  return(results)
}

metrics_huge  <- compute_metrics_extended(hard_huge,  true_graph)
metrics_bd    <- compute_metrics_extended(hard_bd,    true_graph)
metrics_mgm   <- compute_metrics_extended(hard_mgm,   true_graph)
metrics_our1  <- compute_metrics_extended(hard_our1,  true_graph)
metrics_our2  <- compute_metrics_extended(hard_our2,  true_graph)



# Average TPR/FPR/Precision across replicates
summary_metrics_extended <- rbind(
  Huge = colMeans(metrics_huge),
  BDgraph = colMeans(metrics_bd),
  MGM = colMeans(metrics_mgm),
  BMGM_v1 = colMeans(metrics_our1),
  BMGM_v2 = colMeans(metrics_our2)
)
print(round(cbind(summary_metrics_extended, AUC), 3))

sd_metrics_extended <- rbind(
  Huge = apply(metrics_huge,2,sd),
  BDgraph = apply(metrics_bd,2,sd),
  MGM = apply(metrics_mgm,2,sd),
  BMGM_v1 = apply(metrics_our1,2,sd),
  BMGM_v2 = apply(metrics_our2,2,sd)
)
print(round(cbind(sd_metrics_extended, AUC), 3))


#################################################################################

library(BDgraph)
library(dplyr)
library(flux)

evaluate <- function(ppi_df, p, mr_value) {
  l <- choose(p, 2)
  thresholds <- sort(seq(-0.1, 1.1, by = 0.01), decreasing = TRUE)
  nrep <- nrow(ppi_df)

  # Extract matrices
  soft_huge  <- as.matrix(ppi_df[, 1:l])
  soft_bd    <- as.matrix(ppi_df[, (l+1):(2*l)])
  soft_mgm   <- as.matrix(ppi_df[, (2*l+1):(3*l)])
  soft_our1  <- as.matrix(ppi_df[, (3*l+1):(4*l)])
  #soft_our2  <- as.matrix(ppi_df[, (4*l+1):(5*l)])
  true_graph <- as.matrix(ppi_df[, (5*l+1):(6*l)])

  method_list <- list(huge = soft_huge,
                      BDgraph = soft_bd,
                      MGM = soft_mgm,
                      BMGM = soft_our1)
                      #BMGM_v2 = soft_our2)

  roc_data_all <- list()
  auc_summary <- list()

  for (method_name in names(method_list)) {
    pred_matrix <- method_list[[method_name]]
    tpr_mat <- matrix(0, nrow = nrep, ncol = length(thresholds))
    fpr_mat <- matrix(0, nrow = nrep, ncol = length(thresholds))
    auc_vals <- numeric(nrep)

    for (i in 1:nrep) {
      # Convert true edge vector to binary symmetric matrix
      true_mat <- matrix(0, nrow = p, ncol = p)
      true_mat[upper.tri(true_mat)] <- true_graph[i, ]
      true_mat <- true_mat + t(true_mat)
      true_mat <- ifelse(true_mat > 0, 1, 0)

      for (j in seq_along(thresholds)) {
        tr <- thresholds[j]

        # Convert predicted edge vector to symmetric binary matrix
        pred_mat <- matrix(0, nrow = p, ncol = p)
        pred_mat[upper.tri(pred_mat)] <- as.numeric(abs(pred_matrix[i, ]) > tr)
        pred_mat <- pred_mat + t(pred_mat)
        pred_mat <- ifelse(pred_mat > 0, 1, 0)

        res <- BDgraph::compare(pred_mat, true_mat)
        tpr_mat[i, j] <- res["True Positive", 2] / (res["True Positive", 2] + res["False Negative", 2])
        fpr_mat[i, j] <- res["False Positive", 2] / (res["False Positive", 2] + res["True Negative", 2])
      }

      # Sort fpr and tpr for monotonic ROC
      ord <- order(fpr_mat[i, ])
      fpr_mat[i, ] <- fpr_mat[i, ord]
      tpr_mat[i, ] <- tpr_mat[i, ord]

      auc_vals[i] <- flux::auc(fpr_mat[i, ], tpr_mat[i, ])
    }

    # Average and sort again to ensure smooth curve
    mean_fpr <- colMeans(fpr_mat)
    mean_tpr <- colMeans(tpr_mat)
    ord <- order(mean_fpr)

    df <- data.frame(
      fpr = mean_fpr[ord],
      tpr = mean_tpr[ord],
      method = method_name,
      missing_rate = mr_value
    )

    roc_data_all[[method_name]] <- df
    auc_summary[[method_name]] <- c(mean = mean(auc_vals), sd = sd(auc_vals))
  }

  roc_all_df <- dplyr::bind_rows(roc_data_all)
  auc_df <- do.call(rbind, auc_summary)
  auc_df <- data.frame(method = rownames(auc_df), auc_df)
  auc_df$missing_rate <- mr_value

  return(list(roc = roc_all_df, auc = auc_df))
}

p=10
n=500
ppi_df1 <- read.csv("Simulation_bdgraph_p_10_N_500_mr_0.csv")[, -1]
res_00 <- evaluate(ppi_df1, p = 10, mr_value = 0)

ppi_df2 <- read.csv("Simulation_bdgraph_p_10_N_500_mr_0.1.csv")[, -1]
res_01 <- evaluate(ppi_df2, p = 10, mr_value = 0.1)

ppi_df3 <- read.csv("Simulation_bdgraph_p_10_N_500_mr_0.2.csv")[, -1]
res_02 <- evaluate(ppi_df3, p = 10, mr_value = 0.2)

#ppi_df4 <- read.csv("Simulation_bdgraph_p_10_N_500_mr_0.3.csv")[, -1]
#res_03 <- evaluate(ppi_df4, p = 10, mr_value = 0.3)

# Combine for plotting
roc_combined <- bind_rows(res_00$roc, res_01$roc, res_02$roc)#, res_03$roc)

ROC_plot_bd <- ggplot(roc_combined, aes(x = fpr, y = tpr)) +
  geom_line(aes(linetype = method)) +
  scale_linetype_manual(breaks = c("BMGM", "BDgraph", "huge", "MGM"),
                        values=c("solid", "dashed", "dotdash", "dotted")) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color="grey",
              linetype="longdash", size=0.5) +
  labs(x = "False Positive Rate", y = "True Positive Rate", linetype = "") +
  facet_wrap( ~ missing_rate, labeller = as_labeller(appender,
                                                     default = label_parsed))

namefile <- paste0("plot_roc_bdgraph_p_",p,"_n_",n,"_new.png")
ggsave(namefile, ROC_plot_bd, width = 10, height = 3, device = "png", dpi = 600)




evaluation_table <- function(file, p = 10, bfdr = 0.05) {
  obs_graph <- file
  l <- choose(p, 2)

  # Column index matrix for each block
  index <- cbind(seq(1, 5*l, by = l), seq(l, 5*l, by = l))

  # Threshold rules for hard selection
  obs_graph[, index[1,1]:index[1,2]] <- (abs(obs_graph[, index[1,1]:index[1,2]]) >= 0.5) * 1  # Huge
  obs_graph[, index[2,1]:index[2,2]] <- (abs(obs_graph[, index[2,1]:index[2,2]]) >= 0.5) * 1  # BDgraph
  obs_graph[, index[3,1]:index[3,2]] <- (abs(obs_graph[, index[3,1]:index[3,2]]) > 0) * 1     # MGM

  # BMGM: E[FDR] thresholding
  bmgm_ppd <- abs(obs_graph[, index[4,1]:index[4,2]])
  for (i in 1:nrow(bmgm_ppd)) {
    post_inclusion <- bmgm_ppd[i, ]
    fdr_c <- function(c) {
      E_fdr <- sum((1 - post_inclusion) * (post_inclusion > c)) /
        (sum(post_inclusion > c) + rnorm(1, mean = 0, sd = 0.001))
      return(E_fdr)
    }
    pos_c <- seq(0, 1, by = 0.01)
    expected_fdr <- sapply(pos_c, fdr_c)
    pos <- pos_c[expected_fdr < bfdr]
    cutoff <- if (length(pos) > 0) min(pos) else 1
    obs_graph[i, index[4,1]:index[4,2]] <- (bmgm_ppd[i, ] >= cutoff) * 1
  }

  # Initialize metric matrices
  mcc = spec = sens = f1 = tpr = fpr <- matrix(nrow = nrow(obs_graph), ncol = 4)

  for (i in 1:nrow(obs_graph)) {
    adj_template <- diag(p) - diag(p)

    # Convert edge vectors into symmetric matrices
    graphs <- lapply(1:4, function(m) {
      adj <- adj_template
      adj[upper.tri(adj)] <- unlist(obs_graph[i, index[m,1]:index[m,2]])
      adj <- adj + t(adj)
      return(adj)
    })

    # True graph
    true <- adj_template
    true[upper.tri(true)] <- unlist(obs_graph[i, (5*l+1):(6*l)])
    true <- true + t(true)

    # Compute metrics
    comparison <- BDgraph::compare(graphs, true)[,-1]
    mcc[i, ]  <- comparison["MCC", ]
    spec[i, ] <- comparison["Specificity", ]
    sens[i, ] <- comparison["Sensitivity", ]
    f1[i, ]   <- comparison["F1-score", ]
    tpr[i, ]  <- comparison["True Positive", ] / (comparison["True Positive", ] + comparison["False Negative", ])
    fpr[i, ]  <- comparison["False Positive", ] / (comparison["False Positive", ] + comparison["True Negative", ])
  }

  # Combine summaries
  summ <- rbind(
    apply(mcc, 2, mean),
    apply(spec, 2, mean),
    apply(sens, 2, mean),
    apply(f1, 2, mean),
    apply(tpr, 2, mean),
    apply(fpr, 2, mean)
  )

  summ_sd <- rbind(
    apply(mcc, 2, sd),
    apply(spec, 2, sd),
    apply(sens, 2, sd),
    apply(f1, 2, sd),
    apply(tpr, 2, sd),
    apply(fpr, 2, sd)
  )

  colnames(summ) <- colnames(summ_sd) <- c("huge", "BDgraph", "MGM", "BMGM")
  rownames(summ) <- rownames(summ_sd) <- c("MCC", "Specificity", "Sensitivity", "F1", "TPR", "FPR")

  # Format with mean (SD)
  summ[] <- paste0(round(summ, 3), " (", round(summ_sd, 3), ")")

  return(summ)
}

n = 200
mr = 0.3
ppi_df <- read.csv(paste0("Simulation_bdgraph_p_",p,"_N_",n,"_mr_",mr,".csv"))[,-1]
res_00 <- evaluate(ppi_df, p = 10, mr_value = mr)

result_table <- evaluation_table(ppi_df, p = 10)

auc_export <- paste0(round(res_00$auc$mean, 3), " (", round(res_00$auc$sd, 3), ")")

export_tab = t(rbind("AUC" = auc_export,result_table[4,],result_table[1:3,]))
print(export_tab)
