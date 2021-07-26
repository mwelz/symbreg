rm(list = ls()) ; cat("\014")

# load data
load(paste0(getwd(), "/simulations/autocorrelation/sigma05.Rdata"))
results.sigma05 <- measures.overall

load(paste0(getwd(), "/simulations/autocorrelation/sigma2.Rdata"))
results.sigma2 <- measures.overall

data.csv <- rbind(results.sigma05, results.sigma2)
write.csv(data.csv, file = paste0(getwd(), "/simulations/autocorrelation/revised-simulation.csv"))

