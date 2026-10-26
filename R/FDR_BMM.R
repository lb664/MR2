# 135 characters #####################################################################################################################
#' @title False Discovery Rate based on Beta densities finite mixture model
#' @description Post-processing of MR2 Markov chain Monte Carlo output. Estimation of the cutoff of the marginal posterior inclusion
#' probability (mPIP) or posterior edge inclusion probability (PEIP) for a given nominal False Discovery Rate level
#'
#' @param PIP Ordered mPIP or PEIP returned by
#' \code{BMM_EM()}
#' @param w Matrix of allocation probabilities returned by \code{BMM_EM()}
#' @param FDR_nominal Nominal level of False Discovery Rate to be controlled. Default value set at \code{0.05}
#' @param plotting Logical; whether to plot the estimated FDR curve. Default \code{FALSE}
#'
#' @details
#' For details regarding the False Discovery Rate control based on the Beta densities finite mixture model, see
#' \insertCite{Broet2004;textual}{MR2} and \insertCite{Zuber2023;textual}{MR2}
#'
#' @export
#'
#' @return The value returned is a list object \code{list(cutoff, cutoff_idx, FDR_est, FNR_est, opt)}
#' \describe{
#'   \item{\code{cutoff}}{ Cutoff of the mPIP or PEIP for a given nominal FDR level }
#'   \item{\code{cutoff_idx}}{ Index of the ordered mPIP or PEIP that corresponds to the estimated cutoff }
#'   \item{\code{FDR_est}}{ Estimated False Discovery Rate }
#'   \item{\code{FNR_est}}{ Estimated False Negative Rate }
#'   \item{\code{opt}}{ List of options used \code{list(FDR_nominal, plotting)} } }
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
#'
#' PostProc_output <- PostProc(MR2_output, betaHat_Y, betaHat_X)
#'
#' BMM_EM_output <- BMM_EM(PostProc_output$gammaPost)
#' FDR_BMM_output <- FDR_BMM(BMM_EM_output$PIP, BMM_EM_output$w)
#' PostProc_output$thetaPost * (PostProc_output$gammaPost > FDR_BMM_output$cutoff) * 1
#'
#'
#' # Example 2: Analysis of one replication of simulated Scenario IV-Directed pleiotropy with q =
#' # 5 responses, p = 15 exposures and nIV = 100 genetic variants used as IVs and the effect of
#' # the shared pleiotropic pathway on the responses set as \code{2}. The number of expected
#' # exposures directly associated with each response is set at 1 and its variance at 2,
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
#'
#' BMM_EM_output <- BMM_EM(PostProc_output$gammaPost)
#' FDR_BMM_output <- FDR_BMM(BMM_EM_output$PIP, BMM_EM_output$w)
#' PostProc_output$thetaPost * (PostProc_output$gammaPost > FDR_BMM_output$cutoff) * 1
#' }


FDR_BMM <- function(PIP, w, FDR_nominal = 0.05, plotting = FALSE)
{
  N <- length(PIP)

  ########## Local FDR ##########

  ## Reverse order: highest PIP first
  rev_idx <- seq(N, 1)
  w_FPR   <- w[rev_idx, 1]          # P(null | PIP) sorted high to low
  w_FNR   <- 1 - w_FPR

  cum_FDR <- cumsum(w_FPR) / seq_len(N)

  if (sum(cum_FDR < FDR_nominal) > 0)
  {
    cutoff_idx <- max(which(cum_FDR < FDR_nominal))
    cutoff     <- PIP[rev_idx][cutoff_idx]

    FDR_est <- mean(w_FPR[seq_len(cutoff_idx)])
    FNR_est <- mean(w_FNR[(cutoff_idx + 1):N])
  } else {

    if (all(PIP > 0.5))
    {
      cutoff_idx <- 0
      cutoff     <- min(PIP) - 0.01
      FDR_est    <- NA
      FNR_est    <- NA
    } else {
      cutoff_idx <- 1
      cutoff     <- min(PIP) + 0.01
      FDR_est    <- mean(w_FPR[1])
      FNR_est    <- mean(w_FNR[2:N])
    }

  }

  ########## Plotting ##########

  if (isTRUE(plotting) && sum(cum_FDR < FDR_nominal) > 0)
  {
    plot(PIP[rev_idx], cum_FDR,
         xaxt = "n", xlim = c(0, 1), ylim = c(0, 1),
         cex.lab = 1, cex.main = 1,
         xlab = "PIP", ylab = "Nominal and estimated FDR",
         main = "Local and global FDR")
    axis(1, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5))
    abline(h = FDR_nominal, lty = 2, lwd = 2, col = "red")
    abline(v = cutoff,      lty = 2, lwd = 2, col = "red")
  }

  opt <- list(FDR_nominal = FDR_nominal, plotting = plotting)

  output <- list(cutoff     = cutoff,
                 cutoff_idx = cutoff_idx,
                 FDR_est    = FDR_est,
                 FNR_est    = FNR_est,
                 opt        = opt)

  return(output)
}
