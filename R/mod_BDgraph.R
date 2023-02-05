# 135 characters #####################################################################################################################
#' @title Birth and Death sampler of non-decomposable graphs
#' @description Internal function that performs one Markov chain Monte Carlo iteration to sample from the posterior distribution of 
#' a non-decomposable graph. The R code of this function is a modified version of the BDgraph R-package 
#' (\insertCite{Mohammadi2019;textual}{MR2}). For details, see also \insertCite{Mohammadi2015;textual}{MR2} and 
#' \insertCite{Alexopoulos2021;textual}{MR2}
#'
#' @references
#' \insertAllCited{}


########## mod_BDgraph ##########

mod_BDgraph <- function(data, g.start, K.start, samplerPar, 
                        n = NULL, not.cont = NULL, df.prior = 2, g.prior = 0.5, 
                        threshold = 1e-08)
{
  burnin <- samplerPar$BDgraphPar$burnin
  niter <- samplerPar$BDgraphPar$niter
  cores <- samplerPar$BDgraphPar$cores
  
  if (niter < burnin)
  {
    stop("Number of iterations must be larger than the burn-in")
  }
  
  burnin <- floor(burnin)
  cores <- BDgraph::get_cores(cores = cores)
  list_S_n_p <- BDgraph::get_S_n_p(data = data, method = "ggm", n = n, not.cont = not.cont)
  S <- list_S_n_p$S
  n <- list_S_n_p$n
  p <- list_S_n_p$p
  method <- list_S_n_p$method
  colnames_data <- list_S_n_p$colnames_data
  b <- df.prior
  b_star <- b + n
  D <- diag(p)
  Ds <- D + S
  Ts <- chol(solve(Ds))
  Ti <- chol(solve(D))
  g_prior <- BDgraph::get_g_prior(g.prior = g.prior, p = p)
  G <- g.start
  G[g_prior == 1] <- 1
  G[g_prior == 0] <- 0
  G[lower.tri(G, diag(TRUE))] <- 0
  G <- G + t(G)
  K <- K.start
  p_links <- matrix(0, p, p)
  K_hat <- matrix(0, p, p)
  last_graph <- K_hat
  last_K <- K_hat
  
  invisible(capture.output(
  result <- .C("ggm_bdmcmc_ma", as.integer(niter), as.integer(burnin),
               G = as.integer(G), as.double(g_prior), as.double(Ts),
               K = as.double(K), as.integer(p), as.double(threshold),
               K_hat = as.double(K_hat), p_links = as.double(p_links),
               as.integer(b), as.integer(b_star), as.double(Ds), as.integer(niter), PACKAGE = "BDgraph")))
  
  K_hat <- matrix(result$K_hat, p, p, dimnames = list(colnames_data, colnames_data))
  last_graph <- matrix(result$G, p, p, dimnames = list(colnames_data, colnames_data))
  last_K <- matrix(result$K, p, p)
  p_links <- matrix(result$p_links, p, p, dimnames = list(colnames_data, colnames_data))
  p_links[lower.tri(p_links)] <- 0
  output <- list(p_links = p_links, last_K = last_K, last_graph = last_graph)
  
  return(output)
}
