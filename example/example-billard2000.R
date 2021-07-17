rm(list = ls()) ; cat("\014")

# load functions
source(paste0(getwd(), "/funs/symbreg-funs.R"))

# load data from Billard and Diday (2000)
data <- read.csv(paste0(getwd(), "/example/billard2000-data.csv"),
                 header = TRUE, sep = ",", dec = ".")

# Billard and Diday (2000) omit 7th observation
data <- data[-7,] 

# helper objects
X.lower <- cbind(data$x1_lobo, data$x2_lobo)
X.upper <- cbind(data$x1_upbo, data$x2_upbo)

# symbolic regression coefficients
symbolic.regression(X.lower = X.lower, X.upper = X.upper, 
                    Y.lower = data$y_lobo, Y.upper = data$y_upbo)

# boostrap
set.seed(1)
symbolic.regression_bootstrap(X.lower = X.lower, X.upper = X.upper, 
                              Y.lower = data$y_lobo, Y.upper = data$y_upbo, 
                              num.bootstrap.samples = 10000)

# TODO: logical check from Billard and Diday (2000)