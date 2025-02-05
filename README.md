# BMGM: Bayesian Mixed Graphical Model 📊

BMGM (**Bayesian Mixed Graphical Model**) is an R package for **Bayesian inference on mixed graphical models**, allowing the estimation of **conditional dependencies** between **continuous, discrete, categorical, and zero-inflated count** data. The package employs **MCMC sampling** and **spike-and-slab priors** for structure learning and can handle **missing data** during inference.

## 🚀 Installation
To install the package from **GitHub**, run:
```r
install.packages("devtools") # If not installed
devtools::install_github("mauroflorez/BMGM")
```

## Quick Start

```r
library(BMGM)

set.seed(123)
X <- matrix(rnorm(50), 10, 5)
type <- rep("c", 5)  # All continuous variables

fit <- bmgm(X, type, nburn = 500, nsample = 1000)
fit$adj_G
```

## 📖 Main Functions
- `bmgm()`: Fit the Bayesian Mixed Graphical Model  
- `find_lambda()`: Optimize transformation parameter (Arkaprava & Dunson, 2020) 
- `F_transformation()`: Apply transformation for mixed data  

## 📜 License
Released under the **MIT License**. Contributions are welcome!
