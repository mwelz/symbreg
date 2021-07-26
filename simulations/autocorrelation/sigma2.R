rm(list = ls()) ; gc() ; cat("\014")

# load functions
source(paste0(getwd(), "/funs/symbreg-funs.R"))

# scenarios
n.arr       <- 100
alpha.arr   <- 0
beta.arr    <- c(-2, 0, 2) 
rho.arr     <- c(0.2, 0.5)
sigsq.e.arr <- 2
sigsq.z.arr <- 2
sigsq.w.arr <- 2

# settings
R                     <- 1000
num.bootstrap.samples <- 2000
significance.level    <- 0.05

# initialize
set.seed(1832674)
seeds                 <- round(1e7 * runif(R), 0)
significance.quantile <- qnorm(1 - significance.level/2)
count                 <- 1
measures.overall      <- matrix(NA_real_, 6, 7 + 6)
colnames(measures.overall) <- c("n", "alpha", "beta", "rho", 
                                "sigsq_e","sigsq_z", "sigsq_w",
                                "alpha_hat", "SE(alpha_hat)", "alpha_rejection",
                                "beta_hat", "SE(beta_hat)", "beta_rejection")

n       <- 100
alpha   <- 0
beta    <- -2 
rho     <- 0.2
sigsq.e <- 0.5
sigsq.z <- 0.5
sigsq.w <- 0.5
r = 1

for(n in n.arr){
  for(alpha in alpha.arr){
    for(beta in beta.arr){
      for(rho in rho.arr){
        for(sigsq.e in sigsq.e.arr){
          for(sigsq.z in sigsq.z.arr){
            for(sigsq.w in sigsq.w.arr){
              
              # initialize matrix
              measures.across <- matrix(NA_real_, nrow = R, ncol = 6)
              colnames(measures.across) <- c("alpha_hat", "SE(alpha_hat)", "alpha_rejection",
                                             "beta_hat", "SE(beta_hat)", "beta_rejection")
              
              for(r in 1:R){
                
                # set seed
                set.seed(seeds[r])
                
                # sample data
                x    <- rnorm(n, mean = 0, sd = 1)
                eps  <- rnorm(n, mean = 0, sd = sqrt(sigsq.e))
                z.lo <- rnorm(n, mean = 0, sd = sqrt(sigsq.z))
                z.up <- rnorm(n, mean = 0, sd = sqrt(sigsq.z))
                w.lo <- rnorm(n, mean = 0, sd = sqrt(sigsq.w))
                w.up <- rnorm(n, mean = 0, sd = sqrt(sigsq.w))
                
                # DGP
                y <- alpha + beta * x + eps
                
                # generate interval data for y
                y.lo <- y - abs(z.lo)
                y.up <- y + abs(z.up)
                
                # generate interval data for x
                x.lag  <- c(0, x)[1:n]       # x[0] = 0
                x.temp <- x - rho * x.lag 
                x.lo   <- matrix(x.temp - abs(w.lo), ncol = 1)
                x.up   <- matrix(x.temp + abs(w.up), ncol = 1)
                
                # symbolic regression for coefficient estimates
                coef <- symbolic.regression(X.lower = x.lo, 
                                            X.upper = x.up,
                                            Y.lower = y.lo, 
                                            Y.upper = y.up)
                
                # symbolic regression for SEs of coefficient estimates
                se <- symbolic.regression_bootstrap(
                  X.lower = x.lo, 
                  X.upper = x.up,
                  Y.lower = y.lo, 
                  Y.upper = y.up, 
                  num.bootstrap.samples = num.bootstrap.samples)$boostrap.standarderrors
                
                # perform significance test
                significance <- abs(coef / se) > significance.quantile
                
                # store results
                measures.across[r, "alpha_hat"]       <- coef[1]
                measures.across[r, "beta_hat"]        <- coef[2]
                measures.across[r, "SE(alpha_hat)"]   <- se[1]
                measures.across[r, "SE(beta_hat)"]    <- se[2]
                measures.across[r, "alpha_rejection"] <- significance[1]
                measures.across[r, "beta_rejection"]  <- significance[2]
                
              } # FOR r
              
              # store averages across R runs
              measures.overall[count, "n"]       <- n
              measures.overall[count, "alpha"]   <- alpha
              measures.overall[count, "beta"]    <- beta
              measures.overall[count, "rho"]     <- rho
              measures.overall[count, "sigsq_e"] <- sigsq.e
              measures.overall[count, "sigsq_z"] <- sigsq.z
              measures.overall[count, "sigsq_w"] <- sigsq.w
              measures.overall[count, 8:13]      <- colMeans(measures.across)
              
              # update counter
              print(paste0(count, "/6 done at ", Sys.time()))
              count <- count + 1
              
              
            } # FOR sigsq.w
          } # FOR sigsq.z
        } # FOR sigsq.e
      } # FOR rho
    } # FOR beta
  } # FOR alpha
} # FOR n


# store final results
save(measures.overall, file = paste0(getwd(), "/simulations/autocorrelation/sigma2.Rdata"))
