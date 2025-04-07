# 135 characters #####################################################################################################################
#' @title Simulated scenario: Confounding
#' @description Simulator to generate multiple responses and multiple exposures two-sample summary-level data with confounding effects
#'
#' @param nRep Number of simulated replicates
#' @param q Number of simulated summary-level responses
#' @param p Number of simulated summary-level exposures
#' @param nIV Number of genetic variants used as instrumental variables
#' @param nSubject Number of individuals in each sample of individual-level data
#' @param MAF Major allele frequency used in the simulation of the genotypes in each sample individual-level data. Default value set 
#' at \code{0.05}
#' @param thetaRange Sought range of the simulated direct causal effects. Default value set at \code{(-2,2)}
#' @param thetaSparsness Proportion of simulated causal associations between exposures and responses. Default value set at \code{0.30} 
#' @param heritY Heriditability of the responses. Default value set at \code{0.25}
#' @param betaHatRange_X Sought range of the simulated genetic effects on the exposures. Default value set at \code{(-2,2)}
#' @param heritX Heriditability of the exposures. Default value set at \code{0.10}
#' @param rhoX Correlation between the exposures. Default value set at \code{0.60}
#' @param thetaUY Impact of the confounder on the responses. Default value set at \code{1}
#' @param thetaUX Impact of the confounder on the exposures. Default value set at \code{2}
#' @param seed Seed used to initialise the simulation
#'
#' @details
#' For details regarding the simulated scenario, see details in \insertCite{Zuber2023;textual}{MR2}
#'
#' @export
#'
#' @return The value returned is a list object \code{list(betaHat_Y, betaHat_X, theta, par)}
#' \itemize{
#'   \item{\code{betaHat_Y}}{ 3D array (IVs times responses times replicates) of the simulated (summary-level) regression coefficients    
#'         of the  genetic effects on each response }
#'   \item{\code{betaHat_X}}{ 3D array (IVs times exposures times replicates) of the simulated (summary-level) regression coefficients 
#'         of the genetic effects on each exposure }
#'   \item{\code{theta}}{ Matrix ((exposures times replicates) times replicates) of the simulated direct causal effects }
#'   \item{\code{par}}{ List of all parameters used in the simulation \code{list(nRep, q, p, nIV, nSubject, MAF, thetaRange, 
#'         thetaSparsness, heritY, betaHatRange_X, heritX, rhoX, thetaUY, thetaUX, seed) } } }
#'
#' @references
#' \insertAllCited{}

#'
#' @examples
# 100 characters ##################################################################################
#' # Example: Simulation of one replication of simulated Scenario II-Confounding with q = 5 
#' # responses, p = 15 exposures and nIV = 100 genetic variants used as IVs
#'
#' Sim_Confounding <- Scenario_Confounding(nRep = 1, q = 5, p = 15, nIV = 100, seed = 280610971)
#' head(Sim_Confounding$theta)
#' Sim_Confounding$par


Scenario_Confounding <- function(nRep = 2, q = 5, p = 10, nIV = 100, nSubject = 50000, MAF = 0.05, 
                                  thetaRange = c(-2, 2), thetaSparsness = 0.30, heritY = 0.25,  
                                  betaHatRange_X = c(-2, 2), heritX = 0.10, rhoX = 0.60, 
                                  thetaUY = 1, thetaUX = 2, seed = 31122021)
{
  
  if (q == 1) {
    stop("The number of responses should be greater than 1")
  }
  
  ########## Storing ##########
  
  colnames_Y <- paste0("R", seq(1 : q))
  rownames_Y <- paste0("IV", seq(1 : nIV))
  colnames_X <- paste0("E", seq(1 : p))
  rownames_X <- paste0("IV", seq(1 : nIV))
  
  betaHat_YAll <- array(NA, c(nIV, q, nRep), dimnames = list(rownames_Y, colnames_Y, 
                                                             paste0("Rep", seq(1, nRep, 1))))
  betaHat_XAll <- array(NA, c(nIV, p, nRep), dimnames = list(rownames_X, colnames_X, 
                                                             paste0("Rep", seq(1, nRep, 1))))
  thetaAll <- matrix(NA, p * q, nRep, dimnames = list(apply(expand.grid(x = colnames_X, y = colnames_Y), 1, paste, collapse = ","), 
                                                      paste0("Rep", seq(1, nRep, 1))))
  
  ########## Simulating ##########
  
  idxX <- 1 : nSubject
  idxY <- (nSubject + 1) : (2 * nSubject)
  
  for (Rep in 1 : nRep)
  {
    set.seed(seed + Rep)
    cat("Simulating replication:", Rep, "\n")
    
    G <- rbinom(n = nIV * 2 * nSubject, prob = MAF, size = 2)
    G <- matrix(G, nrow = 2 * nSubject, ncol = nIV)
    
    # Confounders
    U <- rnorm(n = 2 * nSubject, mean = 0, sd = 1)
    
    # Create correlated exposures effects
    SigmabetaHat_X <- rhoX ^ (abs(outer(1 : p, 1 : p, "-")))
    SigmaEbetaHat_X <- (betaHatRange_X[2] - betaHatRange_X[1]) / 4 * diag(p)
    SigmaEbetaHat_X <- t(chol(SigmaEbetaHat_X)) %*% SigmabetaHat_X %*% chol(SigmaEbetaHat_X)
    betaHat_X <- mvrnorm(n = nIV, mu = rep(0, times = p), Sigma = SigmaEbetaHat_X)
    
    # Simulate errors exposures
    SigmaEX <- diag(diag(var(G %*% betaHat_X + matrix(thetaUX * U, nrow = 2 * nSubject, ncol = p)))) * (1 - heritX) / heritX
    EX <- mvrnorm(n = 2 * nSubject, mu = rep(0, times = p), Sigma = SigmaEX)
    
    # Simulate exposures
    X <- G %*% betaHat_X + matrix(thetaUX * U, nrow = 2 * nSubject, ncol = p) + EX
    
    # Simulate responses effects
    theta <- mvrnorm(n = p, mu = rep(0, times = q), Sigma = (thetaRange[2] - thetaRange[1]) / 4 * diag(q))
    theta_idx <- matrix(sample(c(rep(1, floor(p * q * thetaSparsness)), 
                        rep(0, times = p * q - floor(p * q * thetaSparsness)))), nrow = p, ncol = q)
    theta <- theta * theta_idx
    thetaAll[, Rep] <- c(theta)
    
    # Simulate responses errors
    rhoY <- 0
    SigmaY <- rhoY ^ (abs(outer(1 : q, 1 : q, "-")))
    SigmaEY <- diag(diag(var(X %*% theta + matrix(thetaUY * U, nrow = 2 * nSubject, ncol = q)))) * (1 - heritY) / heritY
    SigmaEY <- t(chol(SigmaEY)) %*% SigmaY %*% chol(SigmaEY)
    EY <- mvrnorm(n = 2 * nSubject, mu = rep(0, times = q), Sigma = SigmaEY)
    
    # Simulate responses
    Y <- X %*% theta + matrix(thetaUY * U, nrow = 2 * nSubject, ncol = q) + EY
    
    # Create summary-level data
    betaHat_X <- matrix(, nIV, p)
    for (j in 1 : p)
    {
      for(i in 1 : nIV)
      {
        betaHat_X[i, j] <- lm(X[idxX, j] ~ G[idxX, i])$coeff[2]
      }
      cat("Generating betaHat X_j:", j, "\n")
    }
    
    betaHat_Y <- matrix(, nIV, q)
    for (k in 1 : q)
    {
      for(i in 1 : nIV)
      {
        betaHat_Y[i, k] <- lm(Y[idxY, k] ~ G[idxY, i])$coeff[2]
      }
      cat("Generating betaHat Y_k:", k, "\n")
    }
    
    betaHat_XAll[, , Rep] <- betaHat_X
    betaHat_YAll[, , Rep] <- betaHat_Y
  }
  
  if (nRep == 1)
  {
    betaHat_XAll <- betaHat_XAll[, , 1]
    betaHat_YAll <- betaHat_YAll[, , 1]
  }
  
  par = list(nRep = nRep, q = q, p = p, nIV = nIV, nSubject = nSubject, MAF = MAF, 
             thetaRange = thetaRange, thetaSparsness = thetaSparsness, heritY = heritY, 
             betaHatRange_X = betaHatRange_X, heritX = heritX, rhoX = rhoX, 
             thetaUY = thetaUY, thetaUX = thetaUX, seed = seed)
  
  return(list(betaHat_Y = betaHat_YAll, betaHat_X = betaHat_XAll, theta = thetaAll, par = par))
}
