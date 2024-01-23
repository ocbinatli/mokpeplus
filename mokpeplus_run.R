library(pdist)
library(pracma)
library(ManifoldOptim)

source("mokpeplus_projection_train.R")
pdist <- pdist::pdist

args <- commandArgs(trailingOnly = TRUE)
R <- as.numeric(args[[1]])
network <- args[[2]]
tol <- as.numeric(args[[3]])
iter_num <- as.numeric(args[[4]])

nums <- 1:10

data_path <- "./data"
result_path <- "./output"

demo <- function(numb) {
  
  Y <- ...
  K_x <- ...
  K_z <- ...
  
  K_x_feature <- K_x
  K_z_feature <- K_z
  
  N_x <- dim(K_x)[1]
  for(i in 1:N_x) {
    K_x[i, i] <- NA
  }
  
  N_z <- dim(K_z)[1]
  for(i in 1:N_z) {
    K_z[i, i] <- NA
  }
  
  Y[Y == 0] <- NA
  Y <- Y * 0.9  
  
  #set the regularization parameter for cross-domain interactions
  lambda_c <- 1.0
  
  #set the regularization parameters for within-domain similarities
  lambda_x <- 0.1
  lambda_z <- 0.1
  
  #set the kernel width used in the subspace
  sigma_e <- sqrt(R)
  
  #set the parameters list
  parameters <- data.frame(lambda_c, lambda_x, lambda_z, R, numb, sigma_e, network, result_path)
  
  #perform training
  mokpeplus_projection_train(Y, K_x, K_z, K_x_feature, K_z_feature, parameters)
}

parallel::mclapply(nums, demo, mc.cores = 10L)
