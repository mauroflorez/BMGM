
# Application: Kids dataset (Intake and Discharge)

library(mgm)
library(BDgraph)
library(qgraph)
library(readxl)

# Load the dataset
kids <- read_excel("CMC_BMI_12_15_23.xlsx")

# Identify records at intake and discharge
id <- kids$record_id
intake <- which(kids$redcap_event_name %in% c("intake_arm_1", "intake_arm_2"))
discharge <- which(kids$redcap_event_name %in% c("discharge_arm_1", "discharge_arm_2"))

# Create intake dataset
data_intake <- data.frame(
  id = id[intake],
  age = kids$age[intake],
  sex = kids$sex[intake],
  diag = kids$diagnosis[intake],
  edeq_13 = kids$edeq_13[intake],
  edeq_14 = kids$edeq_14[intake],
  edeq_15 = kids$edeq_15[intake],
  edeq_16 = kids$edeq_16[intake],
  edeq_17 = kids$edeq_17[intake],
  edeq_18 = kids$edeq_recode_18[intake],
  edeq_restraint = kids$edeq_restraint[intake],
  edeq_eating = kids$edeq_eating[intake],
  edeq_shape = kids$edeq_shape[intake],
  edeq_weight_con = kids$edeq_weight_con[intake],
  bmi_z = kids$`BMI-Z`[intake]
)

# Keep only anorexia (1) and bulimia (3)
data_intake <- data_intake[data_intake$diag %in% c(1, 3),]
data_intake <- data_intake[!data_intake$id %in% c(319, 271),]  # Remove kids with atypical values
data_intake <- data_intake[rowSums(is.na(data_intake)) < ncol(data_intake),]  # Remove rows with all missing

# Recode diagnosis: anorexia = 0, bulimia = 1
data_intake$diag[data_intake$diag == 1] <- 0
data_intake$diag[data_intake$diag == 3] <- 1

# Create discharge dataset
data_discharge <- data.frame(
  id = id[discharge],
  age = kids$age[discharge],
  sex = kids$sex[discharge],
  diag = kids$diagnosis[discharge],
  edeq_13 = kids$edeq_13[discharge],
  edeq_14 = kids$edeq_14[discharge],
  edeq_15 = kids$edeq_15[discharge],
  edeq_16 = kids$edeq_16[discharge],
  edeq_17 = kids$edeq_17[discharge],
  edeq_18 = kids$edeq_recode_18[discharge],
  edeq_restraint = kids$edeq_restraint[discharge],
  edeq_eating = kids$edeq_eating[discharge],
  edeq_shape = kids$edeq_shape[discharge],
  edeq_weight_con = kids$edeq_weight_con[discharge],
  bmi_z = kids$`BMI-Z`[discharge]
)

# Apply same filtering as intake
data_discharge <- data_discharge[data_discharge$diag %in% c(1, 3),]
data_discharge <- data_discharge[!data_discharge$id %in% c(319),]
data_discharge <- data_discharge[rowSums(is.na(data_discharge)) < ncol(data_discharge),]

data_discharge$diag[data_discharge$diag == 1] <- 0
data_discharge$diag[data_discharge$diag == 3] <- 1

# Final matrices for modeling
data_intake <- as.matrix(data_intake[,-1])
data_discharge <- as.matrix(data_discharge[,-1])
data_discharge <- subset(data_discharge, select = c(1:9, 11, 12, 13, 14, 10))  # match intake structure

# Variable types
type <- c("c", "m", "m", "z", "z", "z", "z", "z", "z", "c", "c", "c", "c", "c")
type2 <- c("g", "c", "c", "p", "p", "p", "p", "p", "p", "g", "g", "g", "g", "g")

# Fit BMGM model to intake
fit1 <- bmgm(data_intake, type, nburn = 5000, nsample = 10000, v_0 = 0.05, v_1 = 1, bfdr = 0.05)
weight <- fit1$adj_Beta
weight[weight < -1] <- -1  # For better visualization

nam <- c("Age", "Sex", "Diag.", "EDE-Q:13", "EDE-Q:14", "EDE-Q:15", "EDE-Q:16", "EDE-Q:17", "EDE-Q:18",
         "EDE-Q:R", "EDE-Q:E", "EDE-Q:S", "EDE-Q:W", "BMI-z")
shape <- c("circle", "triangle", "triangle", rep("square", 6), rep("circle", 5))

colnames(weight) <- nam

# Plot intake graph
qgraph(-weight, labels = nam, shape = shape, theme = "gray", label.cex = 1.3,
       border.width = 2, esize = 5, vsize = 8)

# Fit BMGM model to discharge
fit2 <- bmgm(data_discharge, type, nburn = 5000, nsample = 5000, v_0 = 0.05, v_1 = 1, bfdr = 0.05)
weight2 <- fit2$adj_Beta
weight2[weight2 < -1] <- -1  # For better visualization

# Plot discharge graph
qgraph(-weight2, labels = nam, shape = shape, theme = "gray", label.cex = 1.3,
       border.width = 2, esize = 5, vsize = 8)


