library(pdist)
library(pracma)
library(ManifoldOptim)
source("comp_mokpe_projection_train.R")
pdist <- pdist::pdist

args <- commandArgs(trailingOnly = TRUE)
R <- as.numeric(args[[1]])
network <- args[[2]]

list_vars <- list(f = 1:10, dk = 1:9, tk = 1:9)

demo <- function(list_vars) {
  if (dk == 1) {
    wdd <- sprintf("./kernels/%s_simmat_drugs_aers-bit.txt", network)
  } else if (dk == 2) {
    wdd <- sprintf("./kernels/%s_simmat_drugs_aers-freq.txt", network)
  } else if (dk == 3) {
    wdd <- sprintf("./kernels/%s_simmat_drugs_lambda.txt", network)
  } else if (dk == 4) {
    wdd <- sprintf("./kernels/%s_simmat_drugs_marginalized.txt", network)
  } else if (dk == 5) {
    wdd <- sprintf("./kernels/%s_simmat_drugs_minmaxTanimoto.txt", network)
  } else if (dk == 6) {
    wdd <- sprintf("./kernels/%s_simmat_drugs_sider.txt", network)
  } else if (dk == 7) {
    wdd <- sprintf("./kernels/%s_simmat_drugs_simcomp.txt", network)
  } else if (dk == 8) {
    wdd <- sprintf("./kernels/%s_simmat_drugs_spectrum.txt", network)
  } else if (dk == 9) {
    wdd <- sprintf("./kernels/%s_simmat_drugs_tanimoto.txt", network)
  }
  
  if (tk == 1) {
    wdp <- sprintf("./kernels/%s_simmat_proteins_go.txt", network)
  } else if (tk == 2) {
    wdp <- sprintf("./kernels/%s_simmat_proteins_mismatch-n-k3m1.txt", network)
  } else if (tk == 3) {
    wdp <- sprintf("./kernels/%s_simmat_proteins_mismatch-n-k3m2.txt", network)
  } else if (tk == 4) {
    wdp <- sprintf("./kernels/%s_simmat_proteins_mismatch-n-k4m1.txt", network)
  } else if (tk == 5) {
    wdp <- sprintf("./kernels/%s_simmat_proteins_mismatch-n-k4m2.txt", network)
  } else if (tk == 6) {
    wdp <- sprintf("./kernels/%s_simmat_proteins_ppi.txt", network)
  } else if (tk == 7) {
    wdp <- sprintf("./kernels/%s_simmat_proteins_spectrum-n-k3.txt", network)
  } else if (tk == 8) {
    wdp <- sprintf("./kernels/%s_simmat_proteins_spectrum-n-k4.txt", network)
  } else if (tk == 9) {
    wdp <- sprintf("./kernels/%s_simmat_proteins_sw-n.txt", network)
  }
  set.seed(...)
  
  #within-domain similarity score between drugs
  K_x <- as.matrix(read.table(wdd, header = TRUE))
  K_x <- (K_x + t(K_x)) / 2
  #within-domain similarity score between proteins
  K_z <- as.matrix(read.table(wdp, header = TRUE))
  K_z <- (K_z + t(K_z)) / 2
  
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
  
  fold <- 10
  fold_indices <- 1:N_x
  fold_indices <- mod(fold_indices + fold - 1, fold) + 1
  fold_indices <- fold_indices[randperm(N_x)]
  
  #cross-domain interaction score
  cd <- sprintf("./kernels/%s_admat_dgc.txt", network)
  
  Y <- t(as.matrix(read.table(cd, header = TRUE)))
  Y[Y == 0] <- NA
  Y <- Y * 0.9  
  
  #set the regularization parameter for cross-domain interactions
  lambda_c <- 1.0
  
  #set the regularization parameters for within-domain similarities
  lambda_x <- 0.1
  lambda_z <- 0.1
  
  #set the kernel width used in the subspace
  sigma_e <- sqrt(R)
  
  #set the parameters data frame
  parameters <- data.frame(lambda_c, lambda_x, lambda_z, R, sigma_e, dk, tk, f)
  
  #perform training
  comp_mokpe_projection_train(Y[fold_indices != f,], K_x[fold_indices != f, fold_indices != f], K_z, K_x_feature[fold_indices != f, fold_indices != f], K_z_feature, parameters)

}

parallel::mclapply(list_vars, demo, mc.cores = 52L)