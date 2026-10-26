# 135 characters #####################################################################################################################
#' @title Post-processing of MR2 output
#' @description Post-processing of MR2 Markov chain Monte Carlo output
#'
#' @param output Object returned by \code{MR2()}
#' @param betaHat_Y (number of IVs) times (number of responses) matrix of summary-level responses
#' @param betaHat_X (number of IVs) times (number of exposures) matrix of summary-level exposures
#' @param alpha Area of the tails of a symmetric credible interval (\code{0.05} default) centred around the estimated direct
#' effects and the estimated covariance and (partial) residual correlation matrices
#' @param digits Integer indicating the number of decimal places
#'
#' @details
#' For details regarding the post-processing of the Markov chain Monte Carlo output, see
#' \insertCite{Zuber2023;textual}{MR2}
#'
#' @export
#'
#' @return The value returned is a list object
#' \code{list(gammaPost, thetaPost, thetaPost_CI, DPost, GPost, SigmaPost, SigmaPost_CI,
#' RPost, RPost_CI, Rm1Post, Rm1Post_CI, decomposability_freq, CPO, scaleCPO,
#' freq_outliers, PLML, stability_PLML, DIC, timePostProc, opt)}
#' \describe{
#'   \item{\code{gammaPost}}{ Marginal posterior inclusion probability (mPIP) }
#'   \item{\code{thetaPost}}{ Posterior mean of the direct effects }
#'   \item{\code{thetaPost_CI}}{ Credible interval of the direct effects }
#'   \item{\code{DPost}}{ Posterior mean of the standard deviations of each response }
#'   \item{\code{GPost}}{ Posterior edge inclusion probability (PEIP) }
#'   \item{\code{SigmaPost}}{ Posterior mean of the variance-covariance matrix between the responses after conditioning on the
#'         direct effects }
#'   \item{\code{SigmaPost_CI}}{ (1 - \eqn{\alpha}) credible interval of the variance-covariance matrix }
#'   \item{\code{RPost}}{ Posterior mean of the correlation matrix between the responses }
#'   \item{\code{RPost_CI}}{ (1 - \eqn{\alpha}) credible interval of the correlation matrix }
#'   \item{\code{Rm1Post}}{ Posterior mean of the partial correlation matrix between the responses }
#'   \item{\code{Rm1Post_CI}}{ (1 - \eqn{\alpha}) credible interval of the partial correlation matrix }
#'   \item{\code{decomposability_freq}}{ Frequency of visited decomposable models during the Markov chain Monte Carlo }
#'   \item{\code{CPO}}{ Conditional predictive ordinate for each instrumental variable }
#'   \item{\code{scaleCPO}}{ Scaled conditional predictive ordinate for each instrumental variable }
#'   \item{\code{PLML}}{ Pseudo log-marginal likelihood }
#'   \item{\code{freq_outliers}}{ Frequency of outlier or high-leverage instrumental variables (scaleCPO < 0.01) }
#'   \item{\code{stability_PLML}}{ Stability of the pseudo log-marginal likelihood }
#'   \item{\code{DIC}}{ Deviance Information Criterion }
#'   \item{\code{timePostProc}}{ Time in minutes employed by \code{PostProc} }
#'   \item{\code{opt}}{ List of options used \code{list(alpha, digits)} } }
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \dontrun{
# 100 characters ##################################################################################
#' # Example 1: Analysis of one replication of simulated Scenario II-Confounding with q = 5
#' # responses, p = 15 exposures and nIV = 100 genetic variants used as IVs. The number of
#' # expected exposures directly associated with each response is set at 2 and its variance at 2,
#' # with a priori range of direct association ranging between 0 and 8
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
#' # with a priori range of direct associations ranging between 0 and 7. A non-decomposable
#' # graph for the inverse of the residual correlation between responses is selected
#'
#' Sim_Pleiotropy <- Scenario_Pleiotropy(nRep = 1, q = 5, p = 15, nIV = 100, undirected_A = FALSE,
#'                                       seed = 280610971)
#' betaHat_Y <- Sim_Pleiotropy$betaHat_Y
#' betaHat_X <- Sim_Pleiotropy$betaHat_X
#'
#' MR2_output <- MR2(betaHat_Y, betaHat_X, EVgamma = c(1, 3),
#'                   niter = 15000, burnin = 5000, thin = 10, monitor = 1000,
#'                   nonDecomp = TRUE, seed = 280610971)
#'
#' PostProc_output <- PostProc(MR2_output, betaHat_Y, betaHat_X)
#' }


PostProc <- function(output, betaHat_Y, betaHat_X, alpha = 0.05, digits = 3)
{
  cat("\n")
  cat("##\n")
  cat("## Post-processing of multi-response Mendelian randomisation (MR2)\n")
  cat("## Version: 0.2.0\n")
  cat("##\n")
  cat("## Copyright (C) 2025 L. Bottolo and V. Zuber\n")
  cat("##\n\n")

  startTime <- Sys.time()

  Y <- as.matrix(betaHat_Y)
  X <- as.matrix(betaHat_X)
  rm(betaHat_Y, betaHat_X)

  n <- nrow(Y); m <- ncol(Y); p <- ncol(X)

  if (is.null(colnames(Y))) colnames(Y) <- paste0("R",  seq_len(m))
  if (is.null(rownames(Y))) rownames(Y) <- paste0("IV", seq_len(n))
  if (is.null(colnames(X))) colnames(X) <- paste0("E",  seq_len(p))
  if (is.null(rownames(X))) rownames(X) <- paste0("IV", seq_len(n))

  colnames_Y <- colnames(Y); rownames_Y <- rownames(Y)
  colnames_X <- colnames(X)

  ########## Posterior summaries ##########

  niter <- nrow(output$theta)

  gammaPost <- output$theta
  gammaPost[gammaPost != 0] <- 1
  gammaPost <- .display_round(
    matrix(colMeans(gammaPost), nrow = p, ncol = m,
           dimnames = list(colnames_X, colnames_Y)), digits)

  thetaPost <- .display_round(
    matrix(colMeans(output$theta), nrow = p, ncol = m,
           dimnames = list(colnames_X, colnames_Y)), digits)

  thetaPost_LL <- .display_round(
    matrix(apply(output$theta, 2, quantile, alpha / 2), nrow = p, ncol = m), digits)
  thetaPost_UL <- .display_round(
    matrix(apply(output$theta, 2, quantile, 1 - alpha / 2), nrow = p, ncol = m), digits)
  thetaPost_CI <- array(
    c(thetaPost_LL, thetaPost_UL), dim = c(p, m, 2),
    dimnames = list(colnames_X, colnames_Y, c(alpha / 2, 1 - alpha / 2)))

  DPost <- .display_round(colMeans(output$D), digits)

  GPost <- .display_round(apply(output$G, c(1, 2), mean), digits)
  GPost[upper.tri(GPost, diag = TRUE)] <- NA

  ## ---- Sigma (covariance) posterior ----
  ## Efficient: outer(d, d) * R instead of diag(d) %*% R %*% diag(d)
  SigmaPost_TMP <- matrix(NA_real_, niter, m * m)
  for (iter in seq_len(niter))
  {
    d_iter <- output$D[iter, ]
    SigmaPost_TMP[iter, ] <- c(outer(d_iter, d_iter) * output$R[, , iter])
  }

  SigmaPost <- .display_round(
    matrix(colMeans(SigmaPost_TMP), nrow = m, ncol = m,
           dimnames = list(colnames_Y, colnames_Y)), digits)

  SigmaPost_LL <- .display_round(
    matrix(apply(SigmaPost_TMP, 2, quantile, alpha / 2), nrow = m, ncol = m), digits)
  SigmaPost_LL[upper.tri(SigmaPost_LL)] <- NA

  SigmaPost_UL <- .display_round(
    matrix(apply(SigmaPost_TMP, 2, quantile, 1 - alpha / 2), nrow = m, ncol = m), digits)
  SigmaPost_UL[upper.tri(SigmaPost_UL)] <- NA

  SigmaPost_CI <- array(
    c(SigmaPost_LL, SigmaPost_UL), dim = c(m, m, 2),
    dimnames = list(colnames_Y, colnames_Y, c(alpha / 2, 1 - alpha / 2)))

  ## ---- R (correlation) posterior ----
  RPost <- .display_round(apply(output$R, c(1, 2), mean), digits)
  RPost[upper.tri(RPost, diag = TRUE)] <- NA

  RPost_LL <- .display_round(
    matrix(apply(output$R, c(1, 2), quantile, alpha / 2), nrow = m, ncol = m), digits)
  RPost_LL[upper.tri(RPost_LL, diag = TRUE)] <- NA

  RPost_UL <- .display_round(
    matrix(apply(output$R, c(1, 2), quantile, 1 - alpha / 2), nrow = m, ncol = m), digits)
  RPost_UL[upper.tri(RPost_UL, diag = TRUE)] <- NA

  RPost_CI <- array(
    c(RPost_LL, RPost_UL), dim = c(m, m, 2),
    dimnames = list(colnames_Y, colnames_Y, c(alpha / 2, 1 - alpha / 2)))

  ## ---- Partial correlation (R^{-1}) posterior ----
  valid_iters <- vapply(seq_len(niter), function(iter) {
    all(!is.na(solve(output$R[, , iter])))
  }, logical(1))

  n_valid       <- sum(valid_iters)
  Rm1Post_TMP   <- matrix(NA_real_, n_valid, m * m)
  valid_row     <- 0L
  for (iter in seq_len(niter))
  {
    if (!valid_iters[iter]) next
    valid_row <- valid_row + 1L
    TMP <- -cov2cor(solve(output$R[, , iter]))
    Rm1Post_TMP[valid_row, ] <- c(TMP - 2 * diag(diag(TMP)))
  }

  Rm1Post <- .display_round(
    matrix(colMeans(Rm1Post_TMP), nrow = m, ncol = m,
           dimnames = list(colnames_Y, colnames_Y)), digits)
  Rm1Post[upper.tri(Rm1Post, diag = TRUE)] <- NA

  Rm1Post_LL <- .display_round(
    matrix(apply(Rm1Post_TMP, 2, quantile, alpha / 2), nrow = m, ncol = m), digits)
  Rm1Post_LL[upper.tri(Rm1Post_LL, diag = TRUE)] <- NA

  Rm1Post_UL <- .display_round(
    matrix(apply(Rm1Post_TMP, 2, quantile, 1 - alpha / 2), nrow = m, ncol = m), digits)
  Rm1Post_UL[upper.tri(Rm1Post_UL, diag = TRUE)] <- NA

  Rm1Post_CI <- array(
    c(Rm1Post_LL, Rm1Post_UL), dim = c(m, m, 2),
    dimnames = list(colnames_Y, colnames_Y, c(alpha / 2, 1 - alpha / 2)))

  ########## Checking decomposability ##########

  dec_freq <- mean(vapply(seq_len(niter), function(iter) {
    check_chordal(output$G[, , iter] + diag(m))[[1]]
  }, numeric(1)))

  ########## PLML and DIC ##########

  f_i <- matrix(NA_real_, n, niter)

  for (iter in seq_len(niter))
  {
    theta_iter <- matrix(output$theta[iter, ], nrow = p, ncol = m)
    d_iter     <- output$D[iter, ]
    Sigma_iter <- outer(d_iter, d_iter) * output$R[, , iter]

    for (i in seq_len(n))
      f_i[i, iter] <- mvnfast::dmvn(as.numeric(Y[i, ]),
                                    as.numeric(X[i, ]) %*% theta_iter,
                                    Sigma_iter)
  }

  ICPO     <- rowMeans(1 / f_i)
  CPO      <- 1 / ICPO
  names(CPO) <- rownames_Y
  scaleCPO <- CPO / apply(f_i, 1, max)
  PLML     <- .display_round(sum(log(CPO)), digits)

  freq_outliers   <- .display_round(mean(scaleCPO < 0.01), digits)
  stability_PLML  <- .display_round(mean(log(CPO) < -10), digits)

  bar_D <- mean(-2 * colSums(log(f_i)))

  ## D(theta_bar, Sigma_bar)
  f_i_bar <- vapply(seq_len(n), function(i)
    mvnfast::dmvn(as.numeric(Y[i, ]),
                  as.numeric(X[i, ]) %*% thetaPost,
                  SigmaPost),
    numeric(1))

  D_bar <- -2 * sum(log(f_i_bar))
  DIC   <- .display_round(2 * bar_D - D_bar, digits)

  SigmaPost[upper.tri(SigmaPost)] <- NA

  opt <- list(alpha = alpha, digits = digits)

  endTime     <- Sys.time()
  timePostProc <- paste0(.display_round(as.numeric(difftime(endTime, startTime, units = "min")),
                               digits), "m")

  output <- list(gammaPost          = gammaPost,
                 thetaPost          = thetaPost,
                 thetaPost_CI       = thetaPost_CI,
                 DPost              = DPost,
                 GPost              = GPost,
                 SigmaPost          = SigmaPost,
                 SigmaPost_CI       = SigmaPost_CI,
                 RPost              = RPost,
                 RPost_CI           = RPost_CI,
                 Rm1Post            = Rm1Post,
                 Rm1Post_CI         = Rm1Post_CI,
                 decomposability_freq = dec_freq,
                 CPO                = CPO,
                 scaleCPO           = scaleCPO,
                 freq_outliers      = freq_outliers,
                 PLML               = PLML,
                 stability_PLML     = stability_PLML,
                 DIC                = DIC,
                 timePostProc       = timePostProc,
                 opt                = opt)

  return(output)
}
