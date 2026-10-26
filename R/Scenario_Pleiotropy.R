# 135 characters #####################################################################################################################
#' @title Simulated scenario: Pleiotropy
#' @description Simulator to generate multiple responses and multiple exposures two-sample summary-level data with confounding and
#' pleiotropic effects
#'
#' @param nRep Number of simulated replicates
#' @param q Number of simulated summary-level responses
#' @param p Number of simulated summary-level exposures
#' @param nIV Number of genetic variants used as instrumental variables
#' @param nSubject Number of individuals in each sample of individual-level data
#' @param MAF Minor allele frequency used in the simulation of the genotypes in each sample individual-level data. Default value sets
#' at \code{0.05}
#' @param thetaRange Sought range of the simulated direct effects. Default value sets at \code{(-2,2)}
#' @param thetaSparsness Proportion of simulated direct associations between exposures and responses. Default value sets at \code{0.3}
#' @param herit_Y Heritability of the responses. Default value sets at \code{0.25}
#' @param betaHatRange_X Sought range of the simulated genetic effects on the exposures. Default value sets at \code{(-2,2)}
#' @param herit_X Heritability of the exposures. Default value sets at \code{0.1}
#' @param rho_X Correlation between the exposures. Default value sets at \code{0.6}
#' @param undirected_A Logical parameter for simulating an undirected (\code{TRUE}, default) or directed pleiotropic pathway
#' @param delta_A Effect of the shared pleiotropic pathway on the responses. Default value sets at \code{1}
#' @param psi_Y Impact of the confounder on the responses. Default value sets at \code{1}
#' @param psi_X Impact of the confounder on the exposures. Default value sets at \code{2}
#' @param seed Seed used to initialise the simulation
#'
#' @details
#' For details regarding the simulated scenario, see details in \insertCite{Zuber2023;textual}{MR2}
#'
#' @export
#'
#' @return The value returned is a list object \code{list(betaHat_Y, betaHat_X, theta, parSim)}
#' \describe{
#'   \item{\code{betaHat_Y}}{ 3D array (IVs times responses times replicates) of the simulated (summary-level) regression coefficients
#'         of the  genetic effects on each response }
#'   \item{\code{betaHat_X}}{ 3D array (IVs times exposures times replicates) of the simulated (summary-level) regression coefficients
#'         of the genetic effects on each exposure }
#'   \item{\code{theta}}{ Matrix (exposure-response pairs times replicates) of the simulated direct effects }
#'   \item{\code{parSim}}{ List of all parameters used in the simulation \code{list(nRep, q, p, nIV, nSubject, MAF, thetaRange,
#'         thetaSparsness, herit_Y, betaHatRange_X, herit_X, rho_X, undirected_A, delta_A, psi_Y, psi_X, seed)} } }
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @examples
#' \dontrun{
# 100 characters ##################################################################################
#' # Example 1: Simulation of one replication of simulated Scenario III-Undirected pleiotropy with
#' # q = 5 responses, p = 15 exposures and nIV = 100 genetic variants used as IVs
#'
#' Sim_Pleiotropy <- Scenario_Pleiotropy(nRep = 1, q = 5, p = 15, nIV = 100, seed = 280610971)
#' head(Sim_Pleiotropy$theta)
#' Sim_Pleiotropy$parSim
#'
#'
#' # Example 2: Simulation of one replication of simulated Scenario IV-Directed pleiotropy with
#' # q = 5 responses, p = 15 exposures and nIV = 100 genetic variants used as IVs and the effect
#' # of the shared pleiotropic pathway on the responses set as \code{2}
#'
#' Sim_Pleiotropy <- Scenario_Pleiotropy(nRep = 1, q = 5, p = 15, nIV = 100, undirected_A = FALSE,
#'                                       delta_A = 2, seed = 280610971)
#' head(Sim_Pleiotropy$theta)
#' Sim_Pleiotropy$parSim
#' }


Scenario_Pleiotropy <- function(nRep = 2, q = 5, p = 10, nIV = 100, nSubject = 50000, MAF = 0.05,
                                thetaRange = c(-2, 2), thetaSparsness = 0.3, herit_Y = 0.25,
                                betaHatRange_X = c(-2, 2), herit_X = 0.1, rho_X = 0.6,
                                undirected_A = TRUE, delta_A = 1,
                                psi_Y = 1, psi_X = 2, seed = 31122021)
{

  if (q == 1) {
    stop("The number of responses should be greater than 1")
  }

  ########## Storing ##########

  colnames_Y <- paste0("R", seq_len(q))
  rownames_Y <- paste0("IV", seq_len(nIV))
  colnames_X <- paste0("E", seq_len(p))
  rownames_X <- paste0("IV", seq_len(nIV))

  betaHat_YAll <- array(NA, c(nIV, q, nRep), dimnames = list(rownames_Y, colnames_Y,
                                                             paste0("Rep", seq_len(nRep))))
  betaHat_XAll <- array(NA, c(nIV, p, nRep), dimnames = list(rownames_X, colnames_X,
                                                             paste0("Rep", seq_len(nRep))))
  thetaAll <- matrix(NA, p * q, nRep, dimnames = list(apply(expand.grid(x = colnames_X, y = colnames_Y), 1, paste, collapse = ","),
                                                      paste0("Rep", seq_len(nRep))))

  ########## Simulating ##########

  idxX <- seq_len(nSubject)
  idxY <- (nSubject + 1):(2 * nSubject)

  for (Rep in seq_len(nRep))
  {
    set.seed(seed + Rep)
    cat("Simulating replication:", Rep, "\n")

    G <- rbinom(n = nIV * 2 * nSubject, prob = MAF, size = 2)
    G <- matrix(G, nrow = 2 * nSubject, ncol = nIV)

    # Confounders
    U <- rnorm(n = 2 * nSubject, mean = 0, sd = 1)

    # Create correlated exposures effects
    SigmabetaHat_X <- rho_X ^ (abs(outer(1:p, 1:p, "-")))
    SigmaEbetaHat_X <- (betaHatRange_X[2] - betaHatRange_X[1]) / 4 * diag(p)
    SigmaEbetaHat_X <- t(chol(SigmaEbetaHat_X)) %*% SigmabetaHat_X %*% chol(SigmaEbetaHat_X)
    betaHat_X <- mvrnorm(n = nIV, mu = rep(0, times = p), Sigma = SigmaEbetaHat_X)

    # Simulate errors exposures
    SigmaE_X <- diag(diag(var(G %*% betaHat_X + matrix(psi_X * U, nrow = 2 * nSubject, ncol = p)))) * (1 - herit_X) / herit_X
    E_X <- mvrnorm(n = 2 * nSubject, mu = rep(0, times = p), Sigma = SigmaE_X)

    # Simulate exposures
    X <- G %*% betaHat_X + matrix(psi_X * U, nrow = 2 * nSubject, ncol = p) + E_X

    # Simulate responses effects
    theta <- mvrnorm(n = p, mu = rep(0, times = q), Sigma = (thetaRange[2] - thetaRange[1]) / 4 * diag(q))
    theta_idx <- matrix(sample(c(rep(1, floor(p * q * thetaSparsness)),
                        rep(0, times = p * q - floor(p * q * thetaSparsness)))), nrow = p, ncol = q)
    theta <- theta * theta_idx
    thetaAll[, Rep] <- c(theta)

    # Simulate pleiotropy error
    if (undirected_A == TRUE)
    {
      betaA <- runif(n = nIV, min = betaHatRange_X[1], max = betaHatRange_X[2])
    } else {
      betaA <- runif(n = nIV, min = 0, max = betaHatRange_X[2])
    }
    VA <- var(G %*% betaA) * (1 - herit_X) / herit_X
    EA <- rnorm(n = nSubject, mean = 0, sd = sqrt(VA))

    # Simulate pleiotropic pathway
    A <- G %*% betaA + matrix(rep(EA, 2), ncol = 1, byrow = TRUE)

    # Simulate responses errors
    rho_Y <- 0
    Sigma_Y <- rho_Y ^ (abs(outer(1:q, 1:q, "-")))
    SigmaE_Y <- diag(diag(var(X %*% theta + matrix(psi_Y * U, nrow = 2 * nSubject, ncol = q) + matrix(delta_A * A, 2 * nSubject, q)))) * (1 - herit_Y) / herit_Y
    SigmaE_Y <- t(chol(SigmaE_Y)) %*% Sigma_Y %*% chol(SigmaE_Y)
    E_Y <- mvrnorm(n = 2 * nSubject, mu = rep(0, times = q), Sigma = SigmaE_Y)

    # Simulate responses
    Y <- X %*% theta + matrix(psi_Y * U, nrow = 2 * nSubject, ncol = q) + matrix(delta_A * A, 2 * nSubject, q) + E_Y

    # Create summary-level data (vectorised univariate OLS slopes; see Scenario_Confounding.R for derivation)
    G_X_c <- scale(G[idxX, ], center = TRUE, scale = FALSE)
    X_c   <- scale(X[idxX, ], center = TRUE, scale = FALSE)
    denom_X <- colSums(G_X_c^2)
    betaHat_X <- crossprod(G_X_c, X_c) / denom_X
    cat("Generated betaHat_X for all", p, "exposures\n")

    G_Y_c <- scale(G[idxY, ], center = TRUE, scale = FALSE)
    Y_c   <- scale(Y[idxY, ], center = TRUE, scale = FALSE)
    denom_Y <- colSums(G_Y_c^2)
    betaHat_Y <- crossprod(G_Y_c, Y_c) / denom_Y
    cat("Generated betaHat_Y for all", q, "responses\n")

    betaHat_XAll[, , Rep] <- betaHat_X
    betaHat_YAll[, , Rep] <- betaHat_Y
  }

  if (nRep == 1)
  {
    betaHat_XAll <- betaHat_XAll[, , 1]
    betaHat_YAll <- betaHat_YAll[, , 1]
  }

  parSim <- list(nRep = nRep, q = q, p = p, nIV = nIV, nSubject = nSubject, MAF = MAF,
                 thetaRange = thetaRange, thetaSparsness = thetaSparsness, herit_Y = herit_Y,
                 undirected_A = undirected_A, delta_A = delta_A,
                 betaHatRange_X = betaHatRange_X, herit_X = herit_X, rho_X = rho_X,
                 psi_Y = psi_Y, psi_X = psi_X, seed = seed)

  return(list(betaHat_Y = betaHat_YAll, betaHat_X = betaHat_XAll, theta = thetaAll, parSim = parSim))
}
