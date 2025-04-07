# 135 characters #####################################################################################################################
#' @title Post-processing of MR2 output
#' @description Post-processing of MR2 Markov chain Monte Carlo output
#'
#' @param output Object returned by MR2
#' @param betaHat_Y (number of IVs) times (number of responses) matrix of summary-level responses
#' @param betaHat_X (number of IVs) times (number of exposures) matrix of summary-level exposures
#' @param alpha Area of the tails of a symmetric credible interval (\code{0.05} default) centred around the estimated direct 
#' causal effects and the estimated covariance and (partial) residual correlation matrices
#' @param digits Integer indicating the number of decimal places
#'
#' @details
#' For details regarding the post-processing of the Markov chain Monte Carlo output, see details in 
#' \insertCite{Zuber2023;textual}{MR2}
#'
#' @export
#'
#' @return The value returned is a list object \code{list(gammaPost, thetaPost, thetaPost_CI, DPost, GPost, SigmaPost, SigmaPost_CI, 
#' RPost, RPost_CI, Rm1Post, Rm1Post_CI, decomposability_freq, CPO, scaleCPO, freq_ouliers, PLML, DIC, opt)}
#' \itemize{
#'   \item{\code{gammaPost}}{ Marginal posterior probability of inclusion (mPPI) }
#'   \item{\code{thetaPost}}{ Posterior mean of the direct causal effects }
#'   \item{\code{thetaPost_CI}}{  of the direct causal effects }
#'   \item{\code{DPost}}{ Posterior mean of the standard deviations of each response }
#'   \item{\code{GPost}}{ Edge posterior probability of inclusion (ePPI) }
#'   \item{\code{SigmaPost}}{ Posterior mean of the variance-covariance matrix between the responses after conditioning on the direct 
#'         causal effects }
#'   \item{\code{SigmaPost_CI}}{ (1 - \eqn{\alpha}) credible interval of the variance-covariance matrix between the responses after 
#'         conditioning on the direct causal effects }
#'   \item{\code{RPost}}{ Posterior mean of the correlation matrix between the responses after conditioning on the direct causal 
#'         effects }
#'   \item{\code{RPost_CI}}{ (1 - \eqn{\alpha}) credible interval of the correlation matrix between the responses after conditioning on the 
#'         direct causal effects }
#'   \item{\code{Rm1Post}}{ Posterior mean of the partial correlation matrix between the responses after conditioning on the direct 
#'         causal effects }
#'   \item{\code{Rm1Post_CI}}{ (1 - \eqn{\alpha}) credible interval of the partial correlation matrix between the responses after 
#'         conditioning on the direct causal effects }
#'   \item{\code{decomposability_freq}}{ Frequency of visited decomposability models during the Markov chain Monte Carlo. If 
#'         \code{nonDecomp = TRUE}, frequency can be less than 1 }
#'   \item{\code{CPO}}{ Conditional predictive ordinate for each instrumental variable, see \insertCite{Ntzoufras2008;textual}{MR2} }
#'   \item{\code{scaleCPO}}{ Scaled conditional predictive ordinate for each instrumental variable, see 
#'         \insertCite{Congdon2005;textual}{MR2} }
#'   \item{\code{PLML}}{ Pseudo log-marginal likelihood}
#'   \item{\code{scaleCPO}}{ Scaled conditional predictive ordinate for each instrumental variable, see 
#'         \insertCite{Congdon2005;textual}{MR2} }
#'   \item{\code{freq_ouliers}}{ Frequency outliers or high-leverage and influential instrumental variables defined as the frequency 
#'         of scale CPO < 0.01, see \insertCite{Congdon2005;textual}{MR2} }
#'   \item{\code{stability_PLML}}{ Variance Pseudo log-marginal likelihood }
#'   \item{\code{DIC}}{ Deviation Information Criteria, see \insertCite{Spiegelhalter2002;textual}{MR2} }
#'   \item{\code{timePostProc}}{ Time in minutes employed by \code{PostProc} to perform the analyse }
#'   \item{\code{opt}}{ List of options used \code{list(\eqn{\alpha}, digits)} } }
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
# 100 characters ##################################################################################
#' # Example 1: Analysis of one replication of simulated Scenario II-Confounding with q = 5 
#' # responses, p = 15 exposures and nIV = 100 genetic variants used as IVs. The number of 
#' # expected exposures directly associated with each response is set at 2 and its variance at 2, 
#' # with a priori range of direct causal association ranging between 0 and 8
#'
#' Sim_Confounding <- Scenario_Confounding(nRep = 1, q = 5, p = 15, nIV = 100, seed = 280610971)
#' betaHat_Y <- Sim_Confounding$betaHat_Y
#' betaHat_X <- Sim_Confounding$betaHat_X
#'
#' MR2_output <- MR2(betaHat_Y, betaHat_X, EVgamma = c(2, 2), 
#'                   niter = 7500, burnin = 2500, thin = 5, monitor = 1000, seed = 280610971)
#' PostProc_output <- PostProc(MR2_output, betaHat_Y, betaHat_X)
#'
#'
#' # Example 2: Analysis of one replication of simulated Scenario IV-Directed pleiotropy with q = 
#' # 5 responses, p = 15 exposures and nIV = 100 genetic variants used as IVs and the effect of 
#' # the shared pleiotropic pathway on the responses set as \code{2}. The number of expected 
#' # exposures directly associated with each response is set at 1 and its variance at 3, 
#' # with a priori range of direct causal associations ranging between 0 and 7. A non-decomposable 
#' # graph for the inverse of the residual correlation between responses is selected
#'
#' Sim_Pleiotropy <- Scenario_Pleiotropy(nRep = 1, q = 5, p = 15, nIV = 100, undirectedA = FALSE, 
#'                                       seed = 280610971)
#' betaHat_Y <- Sim_Pleiotropy$betaHat_Y
#' betaHat_X <- Sim_Pleiotropy$betaHat_X
#'
#' MR2_output <- MR2(betaHat_Y, betaHat_X, EVgamma = c(1, 3),
#'                   niter = 15000, burnin = 5000, thin = 10, monitor = 1000, 
#'                   nonDecomp = TRUE, seed = 280610971)
#'
#' PostProc_output <- PostProc(MR2_output, betaHat_Y, betaHat_X)


PostProc <- function(output, betaHat_Y, betaHat_X, alpha = 0.05, digits = 3)
{
  cat("\n")
  cat("##\n")
  cat("## Post-processing of multi-response Medelian randomization (MR2)\n")
  cat("## Version: 0.1.1\n")
  cat("##\n")
  cat("## Copyright (C) 2025 L. Bottolo and V. Zuber\n")
  cat("##\n\n")
  
  start_time <- Sys.time()
  
  Y <- as.matrix(betaHat_Y)
  X <- as.matrix(betaHat_X)
  rm(betaHat_Y)
  rm(betaHat_X)
  
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  p <- dim(X)[2]
  
  if (is.null(colnames(Y)))
  {
    colnames(Y) <- paste0("R", seq(1 : m))
  }
  if (is.null(rownames(Y)))
  {
    rownames(Y) <- paste0("IV", seq(1 : n))
  }
  if (is.null(colnames(X)))
  {
    colnames(X) <- paste0("E", seq(1 : p))
  }
  if (is.null(rownames(X)))
  {
    rownames(X) <- paste0("IV", seq(1 : n))
  }
  
  colnames_Y <- colnames(Y)
  rownames_Y <- rownames(Y)
  colnames_X <- colnames(X)
  rownames_X <- rownames(X)
  
  ########## Posterior summary ##########
  
  niter <- nrow(output$theta)
  thin <- output$samplerPar$thin
  monitor <- output$samplerPar$monitor / thin
  gammaPost <- output$theta
  gammaPost[which(gammaPost != 0)] <- 1
  gammaPost <- apply(gammaPost, 2, mean)
  gammaPost <- round(matrix(gammaPost, nrow = p, ncol = m, dimnames = list(colnames_X, colnames_Y)), digits)
  thetaPost <- apply(output$theta, 2, mean)
  thetaPost <- round(matrix(thetaPost, nrow = p, ncol = m, dimnames = list(colnames_X, colnames_Y)), digits)
  thetaPost_LL <- apply(output$theta, 2, quantile, alpha / 2)
  thetaPost_LL <- round(matrix(thetaPost_LL, nrow = p, ncol = m), digits)
  thetaPost_UL <- apply(output$theta, 2, quantile, 1 - alpha / 2)
  thetaPost_UL <- round(matrix(thetaPost_UL, nrow = p, ncol = m), digits)
  thetaPost_CI <- array(c(thetaPost_LL, thetaPost_UL), dim = c(p, m , 2), dimnames = list(colnames_X, colnames_Y, c(alpha / 2, 1 - alpha / 2)))
  DPost <- round(colMeans(output$D), digits)
  
  GPost <- round(apply(output$G, c(1, 2), mean), digits)
  GPost[upper.tri(GPost, diag = TRUE)] <- NA
  
  SigmaPost_TMP <- NULL
  for (iter in 1 : niter)
  {
    SigmaPost_TMP <- rbind(SigmaPost_TMP, c(diag(output$D[iter, ]) %*% output$R[, , iter] %*% diag(output$D[iter, ])))
  }
  SigmaPost <- round(matrix(colMeans(SigmaPost_TMP), nrow = m, ncol = m, dimnames = list(colnames_Y, colnames_Y)), digits)
  SigmaPost_LL <- apply(SigmaPost_TMP, 2, quantile, alpha / 2)
  SigmaPost_LL <- round(matrix(SigmaPost_LL, nrow = m, ncol = m), digits)
  SigmaPost_LL[upper.tri(SigmaPost_LL, diag = FALSE)] <- NA
  SigmaPost_UL <- apply(SigmaPost_TMP, 2, quantile, 1 - alpha / 2)
  SigmaPost_UL <- round(matrix(SigmaPost_UL, nrow = m, ncol = m), digits)
  SigmaPost_UL[upper.tri(SigmaPost_UL, diag = FALSE)] <- NA
  SigmaPost_CI <- array(c(SigmaPost_LL, SigmaPost_UL), dim = c(m, m , 2), dimnames = list(colnames_Y, colnames_Y, c(alpha / 2, 1 - alpha / 2)))
  
  RPost <- round(apply(output$R, c(1, 2), mean), digits)
  RPost[upper.tri(RPost, diag = TRUE)] <- NA
  RPost_LL <- apply(output$R, c(1, 2), quantile, alpha / 2)
  RPost_LL <- round(matrix(RPost_LL, nrow = m, ncol = m), digits)
  RPost_LL[upper.tri(RPost_LL, diag = TRUE)] <- NA
  RPost_UL <- apply(output$R, c(1, 2), quantile, 1 - alpha / 2)
  RPost_UL <- round(matrix(RPost_UL, nrow = m, ncol = m), digits)
  RPost_UL[upper.tri(RPost_UL, diag = TRUE)] <- NA
  RPost_CI <- array(c(RPost_LL, RPost_UL), dim = c(m, m , 2), dimnames = list(colnames_Y, colnames_Y, c(alpha / 2, 1 - alpha / 2)))
  
  Rm1Post_TMP <- NULL
  for (iter in 1 : niter)
  {
    if (sum(is.na(solve(output$R[, , iter]))) == 0)
    {
      TMP <- -cov2cor(solve(output$R[, , iter]))
      Rm1Post_TMP <- rbind(Rm1Post_TMP, c(TMP - 2 * diag(diag(TMP))))
    }
  }
  Rm1Post <- round(matrix(colMeans(Rm1Post_TMP), nrow = m, ncol = m, dimnames = list(colnames_Y, colnames_Y)), digits)
  Rm1Post[upper.tri(Rm1Post, diag = TRUE)] <- NA
  Rm1Post_LL <- apply(Rm1Post_TMP, 2, quantile, alpha / 2)
  Rm1Post_LL <- round(matrix(Rm1Post_LL, nrow = m, ncol = m), digits)
  Rm1Post_LL[upper.tri(Rm1Post_LL, diag = TRUE)] <- NA
  Rm1Post_UL <- apply(Rm1Post_TMP, 2, quantile, 1 - alpha / 2)
  Rm1Post_UL <- round(matrix(Rm1Post_UL, nrow = m, ncol = m), digits)
  Rm1Post_UL[upper.tri(Rm1Post_UL, diag = TRUE)] <- NA
  Rm1Post_CI <- array(c(Rm1Post_LL, Rm1Post_UL), dim = c(m, m , 2), dimnames = list(colnames_Y, colnames_Y, c(alpha / 2, 1 - alpha / 2)))
  
  ########## Checking decomposability ##########
  
  dec_freq <- 0
  for (iter in 1 : niter)
  {
    dec <- check_chordal((output$G[, , iter] + diag(rep(1, m))))[[1]]
    dec_freq <- dec_freq + dec
    
    # if (iter %% niter == 0)
    # {
    #   cat("Checking decomposability", "\n")
    # }
    
  }
  
  # cat("\n")
  
  dec_freq <- dec_freq / niter
  
  ########## PLML and DIC ##########
  
  f_i <- matrix(NA, n, niter)
  
  for (iter in 1 : niter)
  {
    thetaPost_iter <- matrix(output$theta[iter, ], nrow = p, ncol = m)
    SigmaPost_iter <- diag(output$D[iter, ]) %*% output$R[, , iter] %*% diag(output$D[iter, ])
    
    for (i in 1 : n)
    {
      f_i[i, iter] <- dmvn(as.numeric(Y[i, ]), as.numeric(X[i, ]) %*% thetaPost_iter, SigmaPost_iter)
    }
    
    # if (iter %% niter == 0)
    # {
    #   cat("Computing PLML", "\n")
    # }
    
  }
  
  # cat("\n")
  
  ICPO <- rowMeans(1 / f_i)
  CPO <- 1 / ICPO
  names(CPO) <- rownames_Y
  scaleCPO <- round(CPO / apply(f_i, 1, max), digits)
  PLML <- round(sum(log(CPO)), digits)
  mPLML <- mean(log(CPO))
  stdPLML <- var(log(CPO)) ^(1 /2)
  
  idx <- which(scaleCPO < 0.01)
  freq_ouliers <- round(length(idx) / length(CPO), digits)
  stability_PLML <- round(sum(log(CPO) < -10) / length(CPO), digits)
  
  bar_D <- mean(-2 * colSums(log(f_i)))
  f_i <- rep(NA, n)
  
  for (i in 1 : n)
  {
    f_i[i] <- dmvn(as.numeric(Y[i, ]), as.numeric(X[i, ]) %*% thetaPost, SigmaPost)
  }
  
  D_bar <- -2 * sum(log(f_i))
  pD <- bar_D - D_bar
  DIC <- round(2 * bar_D - D_bar, digits)
  
  SigmaPost[upper.tri(SigmaPost, diag = FALSE)] <- NA
  CPO <- round(CPO, digits)
  scaleCPO <- round(scaleCPO, digits)
  
  opt = list(alpha = alpha, digits = digits)
  
  end_time <- Sys.time()
  timePostProc <- as.numeric(difftime(time1 = end_time, time2 = start_time, units = "min"))
  timePostProc <- paste0(round(timePostProc, digits), "m")
  
  output <- list(gammaPost = gammaPost, 
                 thetaPost = thetaPost, thetaPost_CI = thetaPost_CI, DPost = DPost, 
                 GPost = GPost, SigmaPost = SigmaPost, SigmaPost_CI = SigmaPost_CI, 
                 RPost = RPost, RPost_CI = RPost_CI, Rm1Post = Rm1Post, Rm1Post_CI = Rm1Post_CI, 
                 decomposability_freq = dec_freq, CPO = CPO, scaleCPO = scaleCPO, freq_ouliers = freq_ouliers, 
                 PLML = PLML, stability_PLML = stability_PLML, DIC = DIC, 
                 timePostProc = timePostProc, 
                 opt = opt)
  
  return(output)
}
