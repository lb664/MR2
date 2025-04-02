# 135 characters #####################################################################################################################
#' @title MR2: Multi-response Mendelian randomization
#' @description Markov chain Monte Carlo implementation of Bayesian Multi-response Mendelian randomization model for two-sample 
#' summary-level data
#'
#' @param betaHat_Y (number of IVs) times (number of responses) matrix of summary-level responses
#' @param betaHat_X (number of IVs) times (number of exposures) matrix of summary-level exposures
#' @param EVgamma A priori number of expected exposures directly associated with each response and its variance. For details, see 
#' \insertCite{Kohn2001;textual}{MR2}, \insertCite{Alexopoulos2021;textual}{MR2} and \insertCite{Zuber2023;textual}{MR2}. Default 
#' value is \code{EVgamma = c(2, 2)} which implies a priori range of important exposures associated with each response between 0 and 
#' 6. If a single value is provided between (0, 1), it corresponds to the prior probability of exposure-outcome causal association. 
#' It is compulsory when a single exposure is employed
#' @param tau Prior variance of the (direct causal effect) theta prior with the default value set at \code{1}. Do not change this 
#' value unless prior knowledge of the prior variance is available
#' @param niter Number of Markov chain Monte Carlo iterations (including burn-in)
#' @param burnin Number of Markov chain Monte Carlo iterations to be discarded as burn-in
#' @param thin Store the Markov chain Monte Carlo output at every thin-th iteration
#' @param monitor Display the monitored-th iteration
#' @param fullCov Logical parameter to select full (\code{TRUE}) or sparse (\code{FALSE} default) inverse correlation matrix 
#' Markov chain Monte Carlo sampler
#' @param nonDecomp Logical, \code{TRUE} if a non-decomposable graphical model for the inverse correlation matrix is selected
#' @param alpha Level of significance (\code{0.05} default) to select important exposures for each response by using one at-a-time 
#' response, univariable MR. Used in the construction of the proposal distribution of \code{gamma}
#' @param probAddDel Probability (\code{0.9} default) of adding/deleting one exposure to be directly associated in a model with a 
#' response during Markov chain Monte Carlo exploration
#' @param probAdd Probability (\code{0.9} default) of adding one exposure to be directly associated in a model with a response during 
#' Markov chain Monte Monte Carlo exploration
#' @param probMix Probability (\code{0.25} default) of selecting a geometric proposal distribution that samples the index of the 
#' exposure being added/deleted. For details, see \insertCite{Alexopoulos2021;textual}{MR2}
#' @param geomMean Mean of the geometric proposal distribution that samples the index of the exposure being added/deleted. Default 
#' value is the number of expected exposures directly associated with each response and specified in \code{EVgamma}. For details, 
#' see \insertCite{Alexopoulos2021;textual}{MR2}
#' @param graphParNiter Number of Markov chain Monte Carlo internal iterations to sample the graphical model between responses. 
#' Default value is the max between 16 and the square of the number of responses
#' @param seed Seed used to initialise the Markov chain Monte Carlo algorithm
#'
#' @details
#' For details regarding the model and the algorithm, see \insertCite{Zuber2023;textual}{MR2} and 
#' \insertCite{Alexopoulos2021;textual}{MR2}, respectively
#'
#' @export
#'
#' @return The value returned is a list object \code{list(theta, G, R, D, postMean, hyperPar, samplerPar, opt)}
#' \itemize{
#'   \item{\code{theta}}{ Matrix of the (thinned) samples drawn from the posterior distribution of the direct causal effects }
#'   \item{\code{G}}{ 3D array of the (thinned) samples drawn from the posterior distribution of the graphs }
#'   \item{\code{R}}{ 3D array of the (thinned) samples drawn from the posterior distribution of the correlation matrix between the 
#'         responses }
#'   \item{\code{D}}{ Matrix of the (thinned) samples drawn from the posterior distribution of the standard deviations of each 
#'         responses }
#'   \item{\code{postMean}}{ List of the posterior means of \code{list(gamma, theta, G, R, D)} }
#'   \item{\code{timeMR2}}{ Time in minutes employed by \code{MR2} algorithm to analyse the data }
#'   \item{\code{hyperPar}}{ List of the hyper-parameters \code{list(tau)} and the parameters of the beta-binomial distribution 
#'         derived from \code{EVgamma} }
#'   \item{\code{samplerPar}}{ List of parameters used in the  Markov chain Monte Carlo algorithm 
#'         \code{list(alpha, probAddDel, probAdd, probMix, geomMean, graphParNiter)} }
#'   \item{\code{opt}}{ List of options used \code{list(std, seed)} } }
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
# 100 characters ##################################################################################
#' # Example 1: Analysis of one replication of simulated Scenario II-Confounding with q = 5 
#' # responses, p = 15 exposures and nIV = 100 genetic variants used as IVs. The number of expected 
#' # exposures directly associated with each response is set at 2 and its variance at 2, with a 
#' # priori range of direct causal associations ranging between 0 and 8
#'
#' Sim_Confounding <- Scenario_Confounding(nRep = 1, q = 5, p = 15, nIV = 100, seed = 280610971)
#' betaHat_Y <- Sim_Confounding$betaHat_Y
#' betaHat_X <- Sim_Confounding$betaHat_X
#'
#' MR2_output <- MR2(betaHat_Y, betaHat_X, EVgamma = c(2, 2), 
#'                   niter = 7500, burnin = 2500, thin = 5, monitor = 1000, seed = 280610971)
#' head(MR2_output$postMean$theta)
#'
#'
#' # Example 2: Analysis of one replication of simulated Scenario IV-Directed pleiotropy with q = 
#' # 5 responses, p = 15 exposures and nIV = 100 genetic variants used as IVs and the effect of 
#' # the shared pleiotropic pathway on the responses set as \code{2}. The number of expected 
#' # exposures directly associated with each response is set at 1 and its variance at 2, 
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
#' head(MR2_output$postMean$theta)
#'
#'
#' @importFrom grDevices gray
#' @importFrom graphics abline axis hist lines
#' @importFrom MASS mvrnorm polr
#' @importFrom MCMCpack rwish
#' @importFrom mvnfast dmvn rmvn
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats binomial coefficients cor cov2cor dbeta dbinom dgamma dgeom dnorm glm lm 
#'             kmeans nlm pbinom pgeom pnbinom pnorm quantile qgeom qnorm rbinom rchisq rgamma 
#'             rnorm runif sd var
#' @importFrom truncnorm rtruncnorm
#' @importFrom utils capture.output
#' @importFrom Rdpack reprompt


MR2 <- function(betaHat_Y, betaHat_X, 
                EVgamma = c(2, 2), tau = 1, 
                niter, burnin, thin, monitor, 
                fullCov = FALSE, nonDecomp = FALSE, 
                alpha = 0.05, probAddDel = 0.9, probAdd = 0.5, probMix = 0.25, geomMean = EVgamma[1], 
                graphParNiter = max(16, dim(Y)[2] * dim(Y)[2]), seed = 31122021)
{
  cat("\n")
  cat("##\n")
  cat("## Multi-response Medelian randomization (MR2)\n")
  cat("## Version: 0.1.1\n")
  cat("##\n")
  cat("## Copyright (C) 2025 L. Bottolo and V. Zuber\n")
  cat("##\n")
  cat("## The model selected is:\n")
  cat("##\n")
  
  if (fullCov == FALSE)
  {
    if (nonDecomp)
    {
      cat("## MR2 with non-decomposable graph for unmeasured pleiotropy\n")
    } else {
      cat("## MR2 with decomposable graph for unmeasured pleiotropy\n")
    }
  } else {
    cat("## MR2 with dense graph for unmeasured pleiotropy\n")
  }
  cat("##\n\n")
  
  set.seed(seed)
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
  
  colnames_Y <- colnames(Y)
  rownames_Y <- rownames(Y)
  colnames_X <- colnames(X)
  rownames_X <- rownames(X)
  
  responseType <- rep("Gaussian", m)
  pFixed <- 0
  link <- "probit"
  nCat <- NULL
  cutPoint <- NULL
  negBinomPar <- NULL
  nTrial <- NULL
  
  hyperPar <- NULL
  hyperPar$tau <- tau
  
  samplerPar <- NULL
  samplerPar$niter <- niter
  samplerPar$burnin <- burnin
  samplerPar$thin <- thin
  samplerPar$monitor <- monitor
  samplerPar$fullCov <- fullCov
  samplerPar$nonDecomp <- nonDecomp
  samplerPar$link <- link
  samplerPar$alpha <- alpha
  samplerPar$probAddDel <- probAddDel
  samplerPar$probAdd <- probAdd
  samplerPar$probMix <- probMix
  
  if (length(EVgamma) == 2)
  {
    samplerPar$geomMean <- geomMean
  } else {
    samplerPar$geomMean <- max(1, round(geomMean * p))
  }
  
  samplerPar$GPar$burnin <- 0
  samplerPar$GPar$niter <- graphParNiter
  samplerPar$BDgraphPar$burnin <- 0
  samplerPar$BDgraphPar$niter <- graphParNiter
  samplerPar$BDgraphPar$cores <- 1
  
  ########## Creating matrices to store MCMC output ##########
  
  gammaPropSave <- matrix(NA, niter, p * m, 
                          dimnames = list(paste0("iter", seq(1, niter, 1)), 
                                          apply(expand.grid(x = colnames_X, y = colnames_Y), 1, paste, collapse = ",")))
  BSave <- matrix(NA, (niter - burnin) / thin, p * m, 
                  dimnames = list(paste0("iter", seq(burnin + thin, niter, thin)), 
                                  apply(expand.grid(x = colnames_X, y = colnames_Y), 1, paste, collapse = ",")))
  GSave <- RSave <- array(NA, c(m, m, (niter - burnin) / thin), 
                          dimnames = list(colnames_Y, colnames_Y, 
                                          paste0("iter", seq(burnin + thin, niter, thin))))
  DSave <- matrix(NA, m, (niter - burnin) / thin, 
                  dimnames = list(colnames_Y, 
                                  paste0("iter", seq(burnin + thin, niter, thin))))
  cutPointSave <- list()
  negBinomParSave <- list()
  
  ########## Missing values ##########
  
  if (any(is.na(Y[, which(responseType == "Gaussian")])))
  {
    missC <- TRUE
  } else {
    missC <- FALSE
  }
  
  NA_Y <- NULL
  for (k in 1 : m)
  {
    NA_Y[[k]] <- which(!is.na(Y[, k]))
  }
  
  ########## Computing pvalues needed for gamma proposal ##########
  
  if (pFixed == 0)
  {
    seq_pFixed <- NULL
  } else {
    seq_pFixed <- seq(1, pFixed)
  }
  
  pStar <- p - pFixed
  pval <- matrix(NA, pStar, m)
  
  for (j in 1 : pStar)
  {
    for (k in 1 : m)
    {
      countBinom <- 0
      
      if (responseType[k] == "Gaussian")
      {
        fit <- lm(Y[, k] ~ 0 + X[, c(seq_pFixed, pFixed + j)])
        pval[j, k] <- summary(fit)$coefficients[pFixed + 1, 4]
      } else {
        if (responseType[k] == "binary")
        {
          fit <- glm(factor(Y[, k]) ~ 0 + X[, c(seq_pFixed, pFixed + j)], family = binomial(link = samplerPar$link))
          pval[j, k] <- summary(fit)$coefficients[pFixed + 1, 4]
          Y[, k] <- as.numeric(factor(Y[, k]))
          if (min(Y[, k], na.rm = TRUE) != 0)
          {
            Y[, k] <- Y[, k] - 1
          }
        }
        if (responseType[k] == "ordinal")
        {
          options(warn = -1)
          fit <- polr(factor(Y[, k]) ~ 0 + X[, c(seq_pFixed, pFixed + j)], method = samplerPar$link, Hess = TRUE)
          options(warn = 0)
          coef <- summary(fit)$coefficients
          pval[j, k] <- pnorm(abs(coef[pFixed, 1]), lower.tail = FALSE) * 2
          Y[, k] <- as.numeric(factor(Y[, k]))
          if (min(Y[, k], na.rm = TRUE) != 0)
          {
            Y[, k] <- Y[, k] - 1
          }
        }
        if (responseType[k] == "binomial")
        {
          countBinom <- countBinom + 1
          fit <- glm((Y[, k] / nTrial[countBinom]) ~ 0 + X[, c(seq_pFixed, pFixed + j)], family = binomial(link = samplerPar$link), weights = rep(nTrial[countBinom], n))
          pval[j, k] <- summary(fit)$coefficients[pFixed + 1, 4]
        }
        if (responseType[k] == "negative binomial")
        {
          fit <- glm(Y[, k] ~ 0 + X[, c(seq_pFixed, pFixed + j)], family = "poisson")
          pval[j, k] <- summary(fit)$coefficients[pFixed + 1, 4]
        }
      }
    }
  }
  
  samplerPar$pval <- pval
  samplerPar$mlogpval <- -log10(samplerPar$pval)
  
  ########## Beta-binomial prior gamma ##########
  
  if (length(EVgamma) == 2)
  {  
    Eqg <- EVgamma[1]
    Varqg <- EVgamma[2]
    Epi <- Eqg / pStar
    betaMean <- Epi
    betaVar <- (Varqg - pStar * Epi * (1 - pStar * Epi)) / (pStar * (pStar - 1) * Epi)
    hyperPar$aPi <- (betaMean * betaVar - betaMean) / (betaMean - betaVar)
    hyperPar$bPi <- hyperPar$aPi * (1 - betaMean) / betaMean
    
    if (hyperPar$aPi < 0)
    {
      hyperPar$aPi <- Epi
      hyperPar$bPi <- NA
    }
  
  } else {
    hyperPar$aPi <- EVgamma
    hyperPar$bPi <- NA
  }
  
  mC <- mO <- mOrdinal <- mBinary <- mBinom <- mNegbinom <- 0
  
  for (k in 1 : m)
  {
    if (responseType[k] == "Gaussian")
    {
      mC <- mC + 1
    }
    if (responseType[k] == "binary")
    {
      mBinary <- mBinary + 1
      mO <- mO + 1
    }
    if (responseType[k] == "ordinal")
    {
      mOrdinal <- mOrdinal + 1
      mO <- mO + 1
    }
    if (responseType[k] == "binomial")
    {
      mBinom <- mBinom + 1
      mO <- mO + 1
    }
    if (responseType[k] == "negative binomial")
    {
      mNegbinom <- mNegbinom + 1
      mO <- mO + 1
    }
  }
  
  ########## Initial values gamma ##########
  
  gammaCurr <- rep(0, p * m)
  ind <- rep(NA, m * pFixed)
  if (pFixed > 0)
  {
    for(kk in c(0, 1 : (m - 1)))
    {
      ind[(kk * pFixed + 1) : ((kk + 1) * pFixed)] <- (1 : pFixed) + (kk * pStar + kk * pFixed)
    }
    gammaCurr[ind] <- 1
  }
  for (kk in c(0, 1 : (m - 1)))
  {
    thres <- -log10(samplerPar$alpha / pStar)
    gammaCurr[which(samplerPar$logpval[, kk + 1] > thres) + (kk * pStar + (kk + 1) * pFixed)] <- 1
  }
  
  ########## Initial values decomposable graph, correlation matrix and Beta ##########
  
  GCurr <- diag(rep(1, m))
  diag(GCurr) <- 0
  if ((mBinary > 0) || (mOrdinal > 0))
  {
    cutPointPropLow <- cutPointPropUp <- matrix(NA, n, mBinary + mOrdinal)
  }
  RCurr <- diag(rep(1, m))
  invRCurr <- solve(RCurr)
  BCurr <- rep(0, p * m)
  BCurr[which(gammaCurr != 0)] <- rnorm(length(which(gammaCurr != 0)), 0, 0.1)
  BMat <- matrix(BCurr, m, p, byrow = TRUE)
  d <- rep(1, m)
  D <- diag(d)
  XB <- matrix(NA, n, m)
  for (i in 1 : n)
  {
    for(k in 1 : m)
    {
      XB[i, k] <- sum(X[i, ] * BMat[k, ])
    }
  }
  invSigma <- diag(1 / d) %*% invRCurr %*% diag(1 / d)
  Sigma <- solve(invSigma)
  
  ########## Initial values Z ##########
  
  ZCurr <- matrix(0, n, m)
  if (mC > 0)
  {
    ZCurr[, which(responseType == "Gaussian")] <- Y[, which(responseType == "Gaussian")] - XB[, which(responseType == "Gaussian")]
    ZCurr[, which(responseType == "Gaussian")] <- ZCurr[, which(responseType == "Gaussian")] / sqrt(d[which(responseType == "Gaussian")])
  }
  if (missC == TRUE)
  {
    ZCurr[, which(responseType == "Gaussian")][which(is.na(ZCurr[, which(responseType == "Gaussian")]))] <- rnorm(1)
  }
  
  ########## Initial values cut-points and delta ##########
  
  cutPoint <- list()
  for (k in 1 : m)
  {
    if ((responseType[k] == "Gaussian") | (responseType[k] == "binomial") | (responseType[k] == "negative binomial"))
    {
      cutPoint[[k]] <- NA
    }
    if (responseType[k] == "binary")
    {
      cutPoint[[k]] <- c(-Inf, 0, Inf)
    }
    if (responseType[k] == "ordinal")
    {
      nLev <- max(Y[, k], na.rm = TRUE) + 1
      nCat[k] <- nLev
      cutPoint[[k]] <- rep(NA, nLev + 1)
      cutPoint[[k]][c(1, 2, nLev + 1)] <- c(-Inf, 0, Inf)
      countCat <- 0
      for(kk in 3 : nLev)
      {
        cutPoint[[k]][kk] <- 0.5 + countCat
        countCat <- countCat + 1.5
      }
    }
  }
  WCurr <- ZCurr %*% D
  
  ########## Initial values non-decomposable graph ##########
  
  if (nonDecomp)
  {
    capture.output(
    BD <- mod_BDgraphInit(WCurr, samplerPar))
    invSigma <- BD$last_K
    Sigma <- solve(invSigma)
    d <- sqrt(diag(Sigma))
    RCurr <- diag(1 / d) %*% Sigma %*% diag(1 / d)
    invRCurr <- diag(d) %*% invSigma %*% diag(d)
    D <- diag(d)
    ZCurr <- WCurr %*% diag(1 / d)
  } else {
    GCurr <- diag(m)
  }
  
  ########## Starting MCMC ##########
  
  countIter <- 0
  
  for (iter in 1 : niter)
  {
    Latent <- BVS(Y, NA_Y, X, XB, gammaCurr, BCurr, ZCurr, RCurr, invRCurr, d,
                  responseType, pFixed, 
                  nCat, cutPoint, negBinomPar, 
                  hyperPar, samplerPar)
    
    gammaProp <- Latent[[1]]
    gammaCurr <- Latent[[2]]
    BCurr <- Latent[[3]]
    ZCurr <- Latent[[4]]
    XB <- Latent[[5]]
    cutPoint <- Latent[[6]]
    negBinomPar <- Latent[[7]]
    
    WCurr <- ZCurr %*% D
    
    if (fullCov == FALSE)
    {
      if (nonDecomp)
      {
        capture.output(
        BD <- mod_BDgraph(data = WCurr, g.start = BD$last_graph, K.start = BD$last_K, samplerPar))     
        invSigma <- BD$last_K
        Sigma <- solve(invSigma + diag(rep(10 ^(-10)), m))
        GCurr <- BD$last_graph
      } else {
        GCurr <- Sample_G(WCurr, GCurr, samplerPar)
        HIW <- Sim_HIW(GCurr, t(WCurr) %*% WCurr + diag(rep(1, m)), n + 2)
        Sigma <- HIW[[1]]
        invSigma <- HIW[[2]]
      }
    }
    if (fullCov == TRUE)
    {
      invSigma <- rwish(1 + n + m, solve(t(WCurr) %*% WCurr + diag(rep(1 ,m))))
      Sigma <- solve(invSigma + diag(rep(10 ^(-10), m)))
    }
    d <- sqrt(diag(Sigma))
    RCurr <- diag(1 / d) %*% Sigma %*% diag(1 / d)
    invRCurr <- diag(d) %*% invSigma %*% diag(d)
    D <- diag((d))
    ZCurr <- WCurr %*% diag(1 / d)
    
    gammaPropSave[iter, ] <- gammaProp
    
    if (iter > burnin)
    {
      if (iter %% thin == 0)
      {
        countIter <- countIter + 1
        BSave[countIter, ] <- BCurr
        GSave[, , countIter] <- GCurr
        RSave[, , countIter] <- RCurr
        DSave[, countIter] <- d
        
        cutPointSave[[countIter]] <- cutPoint
        negBinomParSave[[countIter]] <- negBinomPar
      }
    }
    if (iter %% monitor == 0)
    {
      cat(paste("iteration", iter), "\n")
    }
  }
  
  DSave <- t(DSave)
  
  if (all(is.na(unlist(cutPointSave))))
  {
    cutPointSave <- NULL
  }
  
  if (all(is.na(unlist(negBinomParSave))))
  {
    negBinomParSave <- NULL
  }
  
  ########## Gamma proposal parameters ##########
  
  if (is.null(seq_pFixed))
  {
    # samplerPar$gammaProp <- gammaPropSave
    samplerPar$gammaPropMean <- matrix(colMeans(gammaPropSave), p, m)
    colnames(samplerPar$gammaPropMean) <- colnames_Y
    rownames(samplerPar$gammaPropMean) <- colnames_X
    colnames(samplerPar$pval) <- colnames_Y
    rownames(samplerPar$pval) <- colnames_X
    colnames(samplerPar$mlogpval) <- colnames_Y
    rownames(samplerPar$mlogpval) <- colnames_X
  } else {
    # samplerPar$gammaProp <- gammaPropSave
    samplerPar$gammaPropMean <- matrix(colMeans(gammaPropSave), p, m)[-seq_pFixed, ]
    colnames(samplerPar$gammaPropMean) <- colnames_Y
    rownames(samplerPar$gammaPropMean) <- colnames_X[-seq_pFixed]
    colnames(samplerPar$pval) <- colnames_Y
    rownames(samplerPar$pval) <- colnames_X[-seq_pFixed]
    colnames(samplerPar$mlogpval) <- colnames_Y
    rownames(samplerPar$mlogpval) <- colnames_X[-seq_pFixed]
  }
  
  ########## Posterior MCMC ##########
  
  digits <- 3
  postMean <- NULL
  postMean$gamma <- round(matrix(colMeans(abs(BSave) > 0), p, m), digits)
  postMean$theta <- round(matrix(colMeans(BSave), p, m), digits)
  postMean$G <- round(apply(GSave, c(1, 2), mean), digits)
  postMean$R <- round(apply(RSave, c(1, 2), mean), digits)
  postMean$D <- round(colMeans(t(D)), digits)
  
  colnames(postMean$gamma) <- colnames_Y
  rownames(postMean$gamma) <- colnames_X
  colnames(postMean$theta) <- colnames_Y
  rownames(postMean$theta) <- colnames_X
  colnames(postMean$G) <- colnames_Y
  rownames(postMean$G) <- colnames_Y
  colnames(postMean$R) <- colnames_Y
  rownames(postMean$R) <- colnames_Y
  names(postMean$D) <- colnames_Y
  
  opt = list(seed = seed)
  
  end_time <- Sys.time()
  timeMR2 <- as.numeric(difftime(time1 = end_time, time2 = start_time, units = "min"))
  timeMR2 <- paste0(round(timeMR2, digits = 3), "m")
  
  output <- list(theta = BSave, G = GSave, R = RSave, D = DSave, 
                 postMean = postMean, 
                 timeMR2 = timeMR2, 
                 hyperPar = hyperPar,
                 samplerPar = samplerPar, 
                 opt = opt)
  
  return(output)
}


########## BVS ##########

BVS <- function(Y, NA_Y, X, XB, gammaCurr, BCurr, ZCurr, RCurr, invRCurr, d, 
                responseType, pFixed, 
                nCat, cutPoint, negBinomPar, 
                hyperPar, samplerPar)
{
  mlogpval <- samplerPar$mlogpval
  
  aPi <- hyperPar$aPi
  bPi <- hyperPar$bPi
  tau <- hyperPar$tau
  nTrial <- hyperPar$nTrial
  
  gammaProp <- rep(0, length(gammaCurr))
  gammaPropList <- list()
  countBinom <- 0
  countNegBinom <- 0
  
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  p <- dim(X)[2]
  pStar <- p - pFixed
  ind <- 1 : m
  
  mu <- Sigma <- TU <- TL <- matrix(NA, n, m)
  TUProp <- TLProp <- rep(NA, n)
  ZTilde <- ZCurr + XB
  
  for (k in 1 : m)
  {
    NAInd <- NA_Y[[k]] 
    start <- (k - 1) * p + 1
    end <- k * p
    gammaInd <- which(gammaCurr[start : end] > 0)
    
    ########## Proposing gamma ##########
    
    gammaPropList <- Sample_gamma(gammaCurr[(start + pFixed) : end], pStar, mlogpval[, k], samplerPar)
    loggammaPropListRatio <- gammaPropList[[2]]
     
    if (pFixed > 0)
    {
      gammaPropList[[1]] <- c(1 : pFixed, gammaPropList[[1]] + pFixed)
    }
     
    if (length(gammaPropList[[1]] > 0))
    {
      setOne <- gammaPropList[[1]] + (k - 1) * p
      setZero <- ((start : end)[-gammaPropList[[1]]])
    } else {
      setOne <- NULL
      setZero <- (start : end)
    }
    
    if (pFixed > 0)
    {
      gammaStarInit <- gammaPropList[[1]][-c(1 : pFixed)]
    } else {
      gammaStarInit <- gammaPropList[[1]]
    }
    
    if (!is.na(bPi))
    {
      logpriorgammaRatio <- lbeta(length(gammaStarInit) + aPi, pStar - length(gammaStarInit) + bPi) - 
                            lbeta(sum(gammaCurr[start : end]) - pFixed + aPi, pStar - sum(gammaCurr[start : end]) + bPi)
    } else {
      logpriorgammaRatio <- dbinom(length(gammaStarInit), pStar, aPi, log = TRUE) - 
                            dbinom(sum(gammaCurr[start : end]) - pFixed, pStar, aPi, log = TRUE)
    }
    
    if (length(gammaInd) > 0)
    {
      Xgamma <- as.matrix(X[, gammaInd])
      B_k <- BCurr[start : end]
      B_kgamma <- B_k[gammaInd]
      XBgamma <- Xgamma %*% B_kgamma
      tXgamma <- t(Xgamma)
    } else {
      Xgamma <- matrix(0, n, 1)
      B_k <- rep(0, pStar)
      B_kgamma <- 0
      XBgamma <- rep(0, n)
      tXgamma <- t(Xgamma)
    }
    
    if (length(gammaPropList[[1]]) > 0)
    {
      if (pFixed > 0)
      {
        gammaStar <- c(1 : pFixed, gammaStarInit)
      } else {
        gammaStar <- c(gammaStarInit)
      }
      XgammaStar <- as.matrix(X[, gammaStar])
      tXgammaStar <- t(XgammaStar)
    } else {
      if (pFixed > 0)
      {
        gammaStar <- c(1 : pFixed)
        XgammaStar <- as.matrix(X[, gammaStar])
        tXgammaStar <- t(XgammaStar)
      } else {
        gammaStar <- NULL
        XgammaStar <- matrix(0, n, 1)
        tXgammaStar <- t(XgammaStar)
      }
    }
    
    ########## Proposing Beta ##########
    
    RR <- RCurr[k, -k] %*% solve(RCurr[-k, -k])
    mu[, k] <- RR %*% t(ZCurr[, -k])
    Sigma[, k] <- RCurr[k, k] - RR %*% RCurr[-k, k]
    if (responseType[k] == "Gaussian")
    {
      ZTilde[, k] <- d[k] * ZCurr[, k] + XB[, k]
      const <- invRCurr[k, k] / (d[k] * d[k])
      T <- rep(0, n)
      for (kk in ind[-k])
      {
        T <- T + ZCurr[, kk] * invRCurr[k, kk] / d[k]
      }
      if (length(gammaStar) > 1) 
      {
        Q <- diag(rep(1 /(tau), length(gammaStar))) + const * tXgammaStar %*% XgammaStar
        invQ <- solve(Q)
      } else {
        Q <- as.matrix(1 /(tau) ) + const * tXgammaStar %*% XgammaStar
        invQ <- solve(Q)
      }
      muProp <- invQ %*% tXgammaStar %*% (const * ZTilde[, k] + T)
      if (length(gammaInd) > 1)
      {
        Q_OLD <- diag(rep(1 /(tau), length(gammaInd))) + const * tXgamma %*% Xgamma
        invQ_OLD <- solve(Q_OLD)
      } else {
        Q_OLD <- as.matrix(1 /(tau)) + const * tXgamma %*% Xgamma
        invQ_OLD <- solve(as.matrix(1 /(tau)) + const * tXgamma %*% Xgamma)
      }
      muProp_OLD <- invQ_OLD %*% tXgamma %*% (const * ZTilde[, k] + T)
      BProp <- mvrnorm(1, muProp, invQ)
    } else {
      ZTilde[, k] <- ZCurr[, k] + XB[, k]
      const <- invRCurr[k, k]
      T <- rep(0, n)
      for (kk in ind[-k])
      {
        T <- T + ZCurr[, kk] * invRCurr[k, kk]
      }
      if (length(gammaStar) > 1)
      {
        Q <- diag(rep(1 /(tau), length(gammaStar))) + const * tXgammaStar %*% XgammaStar
        invQ <- solve(Q)
      } else {
        Q <- as.matrix(1 /(tau)) + const*tXgammaStar%*%XgammaStar
        invQ <- solve(Q)
      }
      muProp <- invQ %*% tXgammaStar %*% (const * ZTilde[, k] + T)
      if (length(gammaInd) > 1)
      {
        Q_OLD <- diag(rep(1 /(tau), length(gammaInd))) + const * tXgamma %*% Xgamma
        invQ_OLD <- solve(Q_OLD)
      } else {
        Q_OLD <- as.matrix(1 /(tau)) + const * tXgamma %*% Xgamma
        invQ_OLD <- solve(as.matrix(1 /(tau))  + const * tXgamma %*% Xgamma)
      }
      muProp_OLD <- invQ_OLD %*% tXgamma %*% (const * ZTilde[, k] + T)
      BProp <- mvrnorm(1, muProp, invQ)
    }
    if (length(gammaStar) > 0)
    {
      XBProp<-(as.matrix(X[, gammaStar]) %*% BProp)[, 1]
    } else {
      XBProp <- rep(0,n)
    }
    
    ########## Updating Beta and gamma ##########
    
    gammaProp[setOne] <- 1
    gammaProp[setZero] <- 0
    
    if (responseType[k] == "Gaussian")
    {
      likeHolmes <- -0.5 * determinant(Q_OLD, logarithm = TRUE)$modulus + 0.5 * t(muProp_OLD) %*% Q_OLD %*% as.matrix(muProp_OLD) - 0.5 * length(gammaInd) * log(tau)
      likeHolmesProp <- -0.5 * determinant(Q, logarithm = TRUE)$modulus + 0.5 * t(muProp) %*% Q %*% as.matrix(muProp) - 0.5 * length(gammaStar) * log(tau)
      logAccRatio <- likeHolmesProp - likeHolmes + loggammaPropListRatio + logpriorgammaRatio
      if (log(runif(1)) < logAccRatio)
      {
        if (length(gammaStar) > 0)
        {
          BCurr[setOne] <- BProp
          gammaCurr[setOne] <- 1
        }
        gammaCurr[setZero] <- 0
        BCurr[setZero] <- 0
        XBgamma <- XBProp
        XB[, k] <- XBProp
      }
      for (i in 1 : n)
      {
        if (is.na(Y[i, k]))
        {
          ZCurr[i, k] <- rnorm(1, RR %*% (ZCurr[i, -k]), sqrt(RCurr[k, k] - RR %*% RCurr[-k, k]))
        } else {
          ZCurr[i, k] <- (Y[i, k] - XBgamma[i]) / d[k]
        }
      }
    } else {
      if (responseType[k] == "binary")
      {
        likeHolmes <- -0.5 * determinant(Q_OLD, logarithm = TRUE)$modulus + 0.5 * t(muProp_OLD) %*% Q_OLD %*% as.matrix(muProp_OLD) - 0.5*length(gammaInd) * log(tau)
        likeHolmesProp <- -0.5 * determinant(Q, logarithm = TRUE)$modulus + 0.5 * t(muProp) %*% Q %*% as.matrix(muProp) - 0.5 * length(gammaStar) * log(tau)
        logAccRatio <- likeHolmesProp - likeHolmes + loggammaPropListRatio + logpriorgammaRatio
        if (log(runif(1)) < logAccRatio)
        {
          if (length(gammaStar) > 0)
          {
            BCurr[setOne] <- BProp
            gammaCurr[setOne] <- 1
          }
          gammaCurr[setZero] <- 0
          BCurr[setZero] <- 0
          XBgamma <- XBProp
          XB[, k] <- XBProp
        }
        ZCurr[, k] <- ZTilde[, k] - XB[, k]
        TU[, k] <- cutPoint[[k]][Y[, k] + 2] - XBgamma
        TL[, k] <- cutPoint[[k]][Y[, k] + 1] - XBgamma
      }
      
      if (responseType[k] == "ordinal") 
      {
        likeHolmes <- -0.5 * determinant(Q_OLD, logarithm = TRUE)$modulus + 0.5 * t(muProp_OLD) %*% Q_OLD %*% as.matrix(muProp_OLD) - 0.5 * length(gammaInd) * log(tau)
        likeHolmesProp <- -0.5 * determinant(Q, logarithm = TRUE)$modulus + 0.5 * t(muProp) %*% Q %*% as.matrix(muProp) - 0.5 * length(gammaStar) * log(tau)
        logAccRatio <- likeHolmesProp - likeHolmes + loggammaPropListRatio + logpriorgammaRatio
        if (log(runif(1)) < logAccRatio)
        {
          if (length(gammaStar) > 0)
          {
            BCurr[setOne] <- BProp
            gammaCurr[setOne] <- 1
          }
          gammaCurr[setZero] <- 0
          BCurr[setZero] <- 0
          XBgamma <- XBProp
          XB[, k] <- XBProp
        }
        ZCurr[, k] <- ZTilde[, k] - XB[, k]
        TU[, k] <- cutPoint[[k]][Y[, k] + 2] - XBgamma
        TL[, k] <- cutPoint[[k]][Y[, k] + 1] - XBgamma
      }
      
      if (responseType[k] == "binomial")
      {
        countBinom <- countBinom + 1
        TUProp <- qnorm(pbinom(Y[, k], nTrial[countBinom], 1 / (1 + exp(-XBProp))))
        TU[, k] <- qnorm(pbinom(Y[, k], nTrial[countBinom], 1 / (1 + exp(-XBgamma))))
        TLProp <- qnorm(pbinom(Y[, k] - 1, nTrial[countBinom], 1 / (1 + exp(-XBProp))))
        TL[, k] <- qnorm(pbinom(Y[, k] - 1, nTrial[countBinom], 1 / (1 + exp(-XBgamma))))
      }
      
      if (responseType[k] == "negative binomial")
      {
        countNegBinom <- countNegBinom + 1
        TUProp <- qnorm(pnbinom(Y[, k], negBinomPar[countNegBinom], 1 - (1 / (1 + exp(-XBProp)))))
        TU[, k] <- qnorm(pnbinom(Y[, k], negBinomPar[countNegBinom], 1 - (1 / (1 + exp(-XBgamma)))))
        TLProp <- qnorm(pnbinom(Y[, k] - 1, negBinomPar[countNegBinom], 1 - (1 / (1 + exp(-XBProp)))))
        TL[, k]<- qnorm( pnbinom(Y[, k] - 1, negBinomPar[countNegBinom], 1 - (1 / (1 + exp(-XBgamma)))))
        approxInd <- union(which(TU[, k] == TL[, k]), which(TL[, k] == Inf))
        if ((length(approxInd) > 0))
        {
          TL[, k][approxInd] <- qnorm(0.9999999999999)
          TU[, k][approxInd] <- qnorm(0.99999999999999)          
        }
        approxInd <- union(which(TUProp == TLProp), which(TLProp == Inf))
        if ((length(approxInd) > 0))
        {
          TLProp[approxInd] <- qnorm(0.9999999999999)
          TUProp[approxInd] <- qnorm(0.99999999999999)          
        }
      }
      
      if (responseType[k] == "binomial" || responseType[k] == "negative binomial")
      {
        diffLikProp <- pnorm((TUProp[NAInd] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k])) - pnorm((TLProp[NAInd] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k]))
        diffLikProp[which(diffLikProp <= 0)] <- pnorm(Inf) - pnorm(8.2)
        logLikProp <- sum(log(diffLikProp))
        diffLik <- ((pnorm((TU[NAInd, k] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k])) - pnorm((TL[NAInd, k] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k]))))
        diffLik[which(diffLik <= 0)] <- pnorm(Inf) - pnorm(8.2)
        logLik_OLD <- sum(log(diffLik))
        logPrior <- sum(dnorm(BProp, 0, sqrt(tau), log = TRUE)) - sum(dnorm(B_kgamma, 0, sqrt(tau), log = TRUE))
        logAccRatio <- logPrior + dmvnorm(B_kgamma, muProp_OLD, invQ_OLD, log = TRUE) - dmvnorm(BProp, muProp, invQ, log = TRUE) + logpriorgammaRatio
        logAccRatio <- logAccRatio + logLikProp - logLik_OLD + loggammaPropListRatio
        if (log(runif(1)) < logAccRatio)
        {
          if (length(gammaStar) > 0)
          {
            BCurr[setOne] <- BProp
            gammaCurr[setOne] <- 1
            gammaCurr[setZero] <- 0
          }
          BCurr[setZero] <- 0
          XBgamma <- XBProp
          XB[, k] <- XBProp
          TU[, k] <- TUProp
          TL[, k] <- TLProp
        }
      }
      
      ########## Updating cut-points ##########
      
      if (responseType[k] == "ordinal")
      {
        gammaInd <- which(gammaCurr[start : end] > 0)
        if (length(gammaInd) > 0)
        {
          Xgamma <- as.matrix(X[, gammaInd])
          B_k <- BCurr[start : end]
          B_kgamma <- B_k[gammaInd]
          XBgamma <- Xgamma %*% B_kgamma
          tXgamma <- t(Xgamma)
        } else {
          Xgamma <- matrix(0, n, 1)
          B_k <- rep(0, pStar)
          B_kgamma <- 0
          XBgamma <- rep(0, n)
          tXgamma <- t(Xgamma)
        }
        cutPoint_TMP <- rep(NA, length(cutPoint[[k]]) - 2)
        cutPoint_TMP[1] <- 0
        cutPoint_TMP[2 : length(cutPoint_TMP)] <- log(cutPoint[[k]][3 : (length(cutPoint[[k]]) - 1)] - cutPoint[[k]][2 : (length(cutPoint[[k]]) - 2)])
        cutPointProp_TMP <- cutPoint_TMP[2 : length(cutPoint_TMP)] + rnorm(length(2 : length(cutPoint_TMP)), 0, sqrt(0.1))
        cutPointProp <- cutPoint[[k]]
        cutPointProp[3] <- exp(cutPointProp_TMP)[1]
        if (nCat[k] >= 4)
        {
          for (jOrd in 4 : nCat[k])
          {
            cutPointProp[jOrd] <- sum(exp(cutPointProp_TMP)[1 : (jOrd - 2)])
          }
        }
        TUProp <- cutPointProp[Y[, k] + 2] - XBgamma
        TU[, k] <- cutPoint[[k]][Y[, k] + 2] - XBgamma
        TLProp <- cutPointProp[Y[, k] + 1] - XBgamma
        TL[, k] <- cutPoint[[k]][Y[, k] + 1] - XBgamma
        priorRatio <- sum(dnorm(cutPointProp_TMP, 0, sqrt(10), log = TRUE)) - sum(dnorm(cutPoint_TMP[2 : length(cutPoint_TMP)], 0, sqrt(10), log = TRUE))
        diffLikProp <- pnorm((TUProp[NAInd] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k])) - pnorm((TLProp[NAInd] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k]))
        diffLikProp[which(diffLikProp == 0)] < pnorm(Inf) - pnorm(8.2)
        logLikProp <- sum(log(diffLikProp))
        diffLik <- ((pnorm((TU[NAInd, k] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k])) - pnorm((TL[NAInd, k] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k]))))
        diffLik[which(diffLik == 0)] <- pnorm(Inf) - pnorm(8.2)
        logLik_OLD <- sum(log(diffLik))
        likRatio <- logLikProp - logLik_OLD
        
        if (log(runif(1)) < (likRatio + priorRatio))
        {
          cutPoint[[k]][3] <- cutPointProp[3]
          if (nCat[k] >= 4)
          {
            for (jOrd in 4 : nCat[k])
            {
              cutPoint[[k]][jOrd] <- cutPointProp[jOrd]
            }
          }
          TU[, k] <- TUProp
          TL[, k] <- TLProp
        }
      }
      
      ########## Updating negative binomial overdispersion parameter ##########
      
      if (responseType[k] == "negative binomial")
      {
        gammaInd <- which(gammaCurr[start : end] > 0)
        if (length(gammaInd) > 0)
        {
          Xgamma <- as.matrix(X[, gammaInd])
          B_k <-BCurr[start : end]
          B_kgamma <- B_k[gammaInd]
          XBgamma <- Xgamma %*% B_kgamma
          tXgamma <-t(Xgamma)
        } else {
          Xgamma <- matrix(0, n, 1)
          B_k <- rep(0, pStar)
          B_kgamma <- 0
          XBgamma <- rep(0, n)
          tXgamma <- t(Xgamma)
        }
        RR <- RCurr[k, -k] %*% solve(RCurr[-k, -k])
        mu[, k] <- RR %*% t(ZCurr[, -k])
        Sigma[, k] <- RCurr[k, k] - RR %*% RCurr[-k, k]
        rhoProp <- log(negBinomPar[countNegBinom]) + rnorm(1, 0, sqrt(0.05))
        TUProp <- qnorm(pnbinom(Y[, k], exp(rhoProp), 1 - (1 / (1 + exp(-XBgamma)))))
        TLProp <- qnorm(pnbinom(Y[, k] - 1, exp(rhoProp), 1 - (1 / (1 + exp(-XBgamma)))))
        approxInd <- union(which(TUProp == TLProp), which(TLProp == Inf))
        if ((length(approxInd) > 0))
        {
          TLProp[approxInd] <- qnorm(0.9999999999999)
          TUProp[approxInd] <- qnorm(0.99999999999999)
        }
        priorRatio <- dgamma(exp(rhoProp), 2, rate = 1, log = TRUE) - dgamma(negBinomPar[countNegBinom], 2, rate = 1, log = TRUE)
        diffLikProp <- pnorm((TUProp[NAInd] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k])) - pnorm((TLProp[NAInd] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k]))
        diffLikProp[which(diffLikProp == 0)] <- pnorm(Inf) - pnorm(8.2)
        logLikProp <- sum(log(diffLikProp))
        diffLik <- ((pnorm((TU[NAInd, k] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k])) - pnorm((TL[NAInd, k] - mu[NAInd, k]) / sqrt(Sigma[NAInd, k]))))
        diffLik[which(diffLik == 0)] <- pnorm(Inf) - pnorm(8.2)
        logLik_OLD <- sum(log(diffLik))
        likRatio <- logLikProp - logLik_OLD
        if (log(runif(1)) < (likRatio + priorRatio + rhoProp - log(negBinomPar[countNegBinom])))
        {
          negBinomPar[countNegBinom] <- exp(rhoProp)
          TU[, k] <- TUProp
          TL[, k] <- TLProp
        }
      }
      
      ########## Updating Z latent variables ##########
      
      for (i in 1 : n)
      {
        if (!is.na(Y[i, k]))
        {
          ZCurr[i, k] <- rtruncnorm(1, a = TL[i, k], b = TU[i, k], RR %*% (ZCurr[i, -k]), sqrt(RCurr[k, k] - RR %*% RCurr[-k, k]))
          if (is.na(ZCurr[i, k]))
          {
            cat(c(TL[i, k], TU[i, k]), "\n")
          }
        } else {
          ZCurr[i, k] <- rnorm(1, RR %*% (ZCurr[i, -k]), sqrt(RCurr[k, k] - RR %*% RCurr[-k, k]))
        }
      }
    }
  }
  
  output <- list(gammaProp, gammaCurr, BCurr, ZCurr, XB, cutPoint, negBinomPar)
  
  return(output)
}


########## Sample_gamma ##########

Sample_gamma <- function(gammaCurr, p, mlogpval, samplerPar)
{
  probAddDel <- samplerPar$probAddDel
  probAdd <- samplerPar$probAdd
  probMix <- samplerPar$probMix
  geomMean <- samplerPar$geomMean
  
  setOne <- which(gammaCurr != 0)
  setZero <- which(gammaCurr == 0)
  pgam <- length(setOne)
  YXcorr_TMP <- sort(mlogpval, decreasing = TRUE, index.return = TRUE)
  YXcorr_rank_TMP <- YXcorr_TMP$ix
  YXcorr_rank <- intersect(YXcorr_rank_TMP, setOne)
  YXcorr_zero_rank <- intersect(YXcorr_rank_TMP, setZero)
  setOneProp <- setOne
  if ((runif(1) < probAddDel) || (pgam == p) || (pgam == 0))
  {
    if (pgam > 0 & pgam < p)
    {
      if (runif(1) < probAdd)
      {
        moveType <- "add"
        change_j <- Qrnd(p, pgam, probMix, geomMean)
        add_j <- YXcorr_zero_rank[change_j]
        setOneProp <- union(setOne, add_j)
        pgammaProp <- length(setOneProp)
        logPropgamma <- (log(1 - probAdd) - log(probAdd)) + (-log(pgammaProp) - log(Qpdf(change_j, p, pgam, probMix, geomMean)))
      } else {
        moveType <- "delete"
        change_j <- sample.int(pgam, 1)
        del_j <- YXcorr_rank[change_j]
        setOneProp <- setdiff(setOne, del_j)
        pgammaProp <- length(setOneProp)
        logPropgamma <- (log(probAdd) - log(1 - probAdd)) + (log(pgam) + log(Qpdf(change_j, p, pgammaProp, probMix, geomMean)))
      }
    }
    if (pgam == 0 & pgam < p)
    {
      moveType <- "add"
      change_j <- Qrnd(p, 0, probMix, geomMean)
      add_j <- YXcorr_zero_rank[change_j]
      setOneProp <- union(setOne, add_j)
      pgammaProp <- length(setOneProp)
      logPropgamma <- (log(1 - probAdd) - log(1)) + (-log(pgammaProp) - log(Qpdf(change_j, p, pgam, probMix, geomMean)))
    }
    if (pgam > 0 & pgam == p){
      moveType <- "delete"
      change_j <- sample.int(pgam, 1)
      del_j <- YXcorr_rank[change_j]
      setOneProp <- setdiff(setOne, del_j)
      pgammaProp <- length(setOneProp)
      logPropgamma <- (log(probAdd) -log(1) ) + (log(pgam) + log(Qpdf(change_j, p, pgammaProp, probMix, geomMean)))
      }
  } else {
    moveType <- "swap"
    change_j <- Qrnd(p, pgam, probMix, geomMean)
    add_j <- YXcorr_zero_rank[change_j]
    setOneProp <- union(setOne, add_j)
    pgammaProp <- length(setOneProp)
    change_j_OLD <- change_j
    pgammaProp_OLD <- pgammaProp
    cond <- 0
    while (cond == 0)
    {
      change_j <- sample.int(pgammaProp, 1)
      if (abs(change_j - change_j_OLD) > 0)
      {
        cond <- 1
      }
    }
    YX_corr_rank <- intersect(YXcorr_rank_TMP, setOneProp)
    del_j <- YXcorr_rank[change_j]
    setOneProp <- setdiff(setOne, del_j)
    pgammaProp <- length(setOneProp)
    logPropgamma <- (-log(pgammaProp_OLD) - log(Qpdf(change_j_OLD, p, pgam, probMix, geomMean))) + (log(Qpdf(change_j, p, pgammaProp, probMix, geomMean)) + log(pgammaProp_OLD))
  }
  setOneProp <- sort(setOneProp)
  output <- list(setOneProp, logPropgamma, moveType)
  
  return(output)
}

########## Qrnd ##########

Qrnd <- function(p, pgam, u, geomMean)
{
  if (runif(1) < u)
  {
    r <- sample.int(p - pgam, 1)
  } else {
    r <- qgeom(runif(1) * pgeom(p - pgam, 1 /geomMean), 1 /geomMean)
  }
  output <- r
  
  return(output)
}

########## Qpdf ##########

Qpdf <- function(r, p, pgam, u, geomMean)
{
  ff <- u * (1 /(p - pgam)) + (1 - u) * (dgeom(r, 1 /geomMean) / pgeom(p - pgam, 1 /geomMean))
  ouput <- ff
  
  return(ouput)
}


########## Sample_G ##########

Sample_G <- function(W, GCurr, samplerPar, delta = 2)
{
  
  burnin <- samplerPar$GPar$burnin
  niter <- samplerPar$GPar$niter
  
  n <- dim(W)[1]
  J <- dim(W)[2]
  Phi <- diag(rep(1, J))
  WW <- t(W) %*% W
  
  G <- generate_graph_zo(GCurr, delta, Phi, n, WW, burnin, niter)[[1]]
  output <- G
  
  return(output)
}

########## ripcliques_to_jtree_zo ##########

ripcliques_to_jtree_zo <- function(cliques)
{
  p <- dim(cliques)[1]
  clique_sizes <- rep(0, p)
  clique_sizes <- colSums(cliques)
  num_cliques <- length(which(clique_sizes != 0))
  score <- rep(0, p)
  jtree <- matrix(0, p, p)
  if (num_cliques >= 2)
  {
    for (i in 2 : num_cliques)
    {
      for (k in 1 : (i - 1))
      {
        score[k] <- sum(intersect_zo(cliques[, i], cliques[, k]))
      }
      if (max(score) != 0)
      {
        argmax <- which(score == max(score))[1]
        j <- argmax
        jtree[i, j] <- 1
        jtree[j, i] <- 1
      }
    }
  }
  output <- jtree
  
  return(output)
}

########## sepsize_zo ##########

sepsize_zo <- function(cliques, jtree)
{
  p <- dim(cliques)[1]
  sepsize <- matrix(0, p, p)
  clique_sizes <- rep(0, p)
  clique_sizes <- colSums(cliques)
  num_cliques <- length(which(clique_sizes != 0))
  if (num_cliques >= 1)
  {
    for (i in 1 : num_cliques)
    {
      if ((i + 1) <= num_cliques)
      {
        for (j in (i + 1) : num_cliques)
        {
          if (jtree[i,j] == 1)
          {
            sepsize[i, j] <- sum(intersect_zo(cliques[, i], cliques[, j]))
            sepsize[j, i] <- sepsize[i, j]
          }
        }
      }
    }
  }
  output <- sepsize
  
  return(output)
}

########## reachability_graph ##########

reachability_graph <- function(g)
{
  p <- dim(g)[1]
  A <- g
  reach_graph <- matrix(0, p, p)
  for (i in 1 : (p - 1))
  {
    reach_graph <- reach_graph + A
    A <- A %*% g
  }
  reach_graph[which(reach_graph > 0)] <- 1
  output <- reach_graph
  
  return(output)
}

########## next_edge_candidate ##########

next_edge_candidate <- function(g_current)
{
  p <- max(dim(g_current)[1], dim(g_current)[2])
  iota <- sample(p)
  i <- iota[1]
  j <- iota[2]
  edge_candidate <- c(i, j)
  output <- edge_candidate
  
  return(output)
}

########## check_edge_delete_zo ##########

check_edge_delete_zo <- function(i, j, cliques)
{
  p <- dim(cliques)[1]
  delete_ok <- 1
  C <- rep(0, p)
  count_cliques <- 0
  clique_sizes <- rep(0, p)
  clique_sizes <- colSums(cliques)
  num_cliques <- length(which(clique_sizes != 0))
  edge_col_vec <- rep(0, p)
  edge_col_vec[i] <- 1
  edge_col_vec[j] <- 1
  if (num_cliques >= 1)
  {
    for (k in 1 : num_cliques)
    {
      clique_k <- cliques[, k]
      if (sum(intersect_zo(edge_col_vec, clique_k)) == 2)
      {
        C <- clique_k
        count_cliques <- count_cliques + 1
        if (count_cliques > 1) 
        {
          delete_ok <- 0
          C <- rep(9, 3)
          
          break
        }
      }
    }
  }
  output <- list(delete_ok, C)
  
  return(output)
}

########## check_edge_delete_zo ##########

find_clique_containing_zo<-function(a, cliques){

  p <- dim(cliques)[1]
  all_indicies <- rep(0, p)
  all_indicies <- which(cliques[a, ]!=0)
  index_a <- min(all_indicies)
  output <- index_a
  
  return(output)
}

########## neighbours_vertex_zo ##########

neighbours_vertex_zo <- function(g, v_i)
{
  p <- dim(g)[1]
  ns <- rep(0, p)
  ns <- g[, v_i]
  output <- ns
  
  return(output)
}

########## parents_vertex_zo ##########

parents_vertex_zo <- function(g, i)
{
  p <- dim(g)[1]
  ps <- rep(0,p)
  nbs <- neighbours_vertex_zo(g, i);
  ps <- nbs
  ps[i : p] <- 0
  output <- ps
  
  return(output)
}

########## check_edge_add_same_component_zo ##########

check_edge_add_same_component_zo <- function(a, b, jtree, sepsize, cliques)
{
  edge_ok <- 0
  quit <- 0
  fork <- 0
  CASE_same_branch <- 0
  CASE_sats_on_b2_branch <- 0
  locate_a2 <- 0
  p <- dim(cliques)[1]
  ab_edge_vec <- rep(0, p)
  ab_edge_vec[a] <- 1
  ab_edge_vec[b] <- 1
  index_b <- find_clique_containing_zo(b, cliques)
  index_a <- find_clique_containing_zo(a, cliques)
  index_b2 <- max(c(index_a, index_b))
  if (index_b2 == index_b)
  {
    b2 <- b
    a2 <- a
  }else{
    b2 <- a
    a2 <- b
  }
  index <- index_b2
  next_index <- which(parents_vertex_zo(jtree, index) !=0 )
  if ((index_b2 - 1) >= 1)
  {
    for (dummy in 1 : (index_b2 - 1))
    {
      if (sum(next_index) == 0 & locate_a2 == 0)
      {
        CASE_same_branch <- 0
        
        break
      }
      clique_next_index <- (cliques[, next_index])
      if (clique_next_index[a2] == 1 & sum(next_index) != 0)
      {
        index_a2 <- next_index
        locate_a2 <- 1
        CASE_same_branch <- 1
        
        break
      }
      index <- next_index
      next_index <- which(parents_vertex_zo(jtree, index) != 0)
    }
  }
  if (CASE_same_branch == 0)
  {
    index_a2 <- find_clique_containing_zo(a2, cliques)
  }
  clique_a2_int_clique_b2 <- intersect_zo(cliques[, index_a2], cliques[, index_b2])
  s <- sum(clique_a2_int_clique_b2)
  if (s == 0)
  {
    edge_ok <- 0
    C <- rep(0, p)
    quit <- 1
  }
  
  if (jtree[index_a2, index_b2] == 1)
  {
    edge_ok <- 1
    C <- union_zo(ab_edge_vec, clique_a2_int_clique_b2)
    quit <- 1
  }
  if (quit == 0)
  {
    if (CASE_same_branch == 1)
    {
      bottom <- index_b2
      next_parent_b2 <- which(parents_vertex_zo(jtree, index_b2) != 0)
      while ((sum(next_parent_b2) > 0) & (next_parent_b2 >= index_a2))
      {
        if (sepsize[next_parent_b2, bottom] == s)
        {
          edge_ok <- 1
          C <- union_zo(ab_edge_vec, clique_a2_int_clique_b2)
          quit <- 1
          
          break
        }
        bottom <- next_parent_b2
        next_parent_b2 <- which(parents_vertex_zo(jtree, bottom) !=0 )
        if (length(next_parent_b2) == 0)
        {
          
          break
        }
      }
    }else{
      next_parent_a2 <- which(parents_vertex_zo(jtree, index_a2) != 0)
      next_parent_b2 <- which(parents_vertex_zo(jtree, index_b2) != 0)
      ancestors_a2 <- rep(0, index_a2)
      ancestors_b2 <- rep(0, index_b2)
      ancestors_a2_col_vec <- rep(0, p)
      ancestors_b2_col_vec <- rep(0, p)
      next_ancestor_a2 <- index_a2
      next_ancestor_b2 <- index_b2
      if (index_a2 >= 1)
      {
        for (count in 1 : index_a2)
        {
          ancestors_a2[count] <- next_ancestor_a2
          ancestors_a2_col_vec[next_ancestor_a2] <- 1
          next_ancestor_a2 <- which(parents_vertex_zo(jtree, next_ancestor_a2) != 0)
          if (sum(next_ancestor_a2) == 0)
          {
            
            break
          }
        }
      }
      if (index_b2 >= 1)
      {
        for (count in 1 : index_b2)
        {
          ancestors_b2[count] <- next_ancestor_b2
          ancestors_b2_col_vec[next_ancestor_b2] <- 1
          next_ancestor_b2 <- which(parents_vertex_zo(jtree, next_ancestor_b2) != 0)
          if (sum(next_ancestor_b2) == 0)
          {
            
            break
          }
        }
      }
      fork_set <- intersect_zo(ancestors_a2_col_vec, ancestors_b2_col_vec)
      fork <- max(which((fork_set != 0)))
      num_seps_branch_b2 <- length(which(ancestors_b2 > fork))
      if (num_seps_branch_b2 > 0)
      {
        for (count in 1 : num_seps_branch_b2)
        {
          index_row <- ancestors_b2[count]
          index_col <- ancestors_b2[count + 1]
          if (sepsize[index_row, index_col] == s)
          {
            edge_ok <- 1
            C <- union_zo(ab_edge_vec, clique_a2_int_clique_b2)
            quit <- 1
            CASE_sats_on_b2_branch <- 1
            
            break
          }
        }
      }
      if (quit == 0 & CASE_sats_on_b2_branch == 0)
      {
        num_seps_branch_a2 <- length(which(ancestors_a2 > fork))
        if (num_seps_branch_a2 > 0)
        {
          for (count in 1 : num_seps_branch_a2)
          {
            index_row <- ancestors_a2[count + 1]
            index_col <- ancestors_a2[count]
            if (sepsize[index_row, index_col] == s)
            {
              edge_ok <- 1
              C <- union_zo(ab_edge_vec, clique_a2_int_clique_b2)
              quit <- 1
              
              break
            }
          }
        }
      }
    }
  }
  if (quit == 0)
  {
    edge_ok <- 0
    C <- NULL
  }
  output <- list(edge_ok, C)
  
  return(output)
}

########## next_graph_candidate_zo ##########

next_graph_candidate_zo <- function(g_current, jtree, sepsize, cliques, reach_graph)
{
  p <- max(dim(g_current)[1], dim(g_current)[2])
  g_proposal <- g_current
  CASE_add <- 0
  CASE_delete <- 0
  while(identical(g_proposal, g_current))
  {
    edge <- next_edge_candidate(g_current)
    i <- edge[1]
    j <- edge[2]
    if (g_current[i, j] == 1)
    {
      ss <- check_edge_delete_zo(i, j, cliques)
      CASE_delete <- ss[[1]]
      C_potential_delete <- ss[[2]]
      if (CASE_delete == 1)
      {
        g_proposal[i, j] <- 0
        g_proposal[j, i] <- 0
        CASE_delete <- 1
        a <- i
        b <- j
        C <- 
        C_potential_delete
      }
    }else{
      if (reach_graph[i, j] == 0)
      {
        CASE_add <- 1
        a <- i
        b <- j
        C <- c(a, b)
        g_proposal[i, j] <- 1
        g_proposal[j, i] <- 1
      }else{
        ss <- check_edge_add_same_component_zo(i, j, jtree, sepsize, cliques)
        CASE_add <- ss[[1]]
        C_potential_add <- ss[[2]]
        if (CASE_add == 1)
        {
          a <-i
          b <- j
          C <- C_potential_add
          g_proposal[i, j] <- 1
          g_proposal[j, i] <- 1
        }
      }
    }
  }
  output <- list(g_proposal, a, b, C, CASE_add, CASE_delete)
  
  return(output)
}

########## ln_h_ratio_zo ##########

ln_h_ratio_zo <- function(C, a, b, delta, Phi)
{
  C_nodes <- which(C != 0)
  D <- c(a, b)
  Sq2 <- setdiff(C_nodes, D)
  numSq2 <- length(Sq2)
  delta_star <- delta - 1
  Phi_CC <- matrix(0, length(C_nodes), length(C_nodes))
  Phi_CC <- rbind(cbind(matrix(Phi[Sq2, Sq2], length(Sq2), length(Sq2)), matrix(Phi[Sq2,D], length(Sq2), length(D))), 
                  cbind(matrix(Phi[D, Sq2], length(D), length(Sq2)), matrix(Phi[D, D], length(D), length(D))))
  L <- t(chol(Phi_CC))
  L_DD <- L[c(numSq2 + 1, numSq2 + 2), c(numSq2 + 1, numSq2 + 2)]
  l_aa <- L_DD[1, 1]
  l_bb <- L_DD[2, 2]
  l_ab <- L_DD[2, 1]
  ln_top <- (delta_star + numSq2 + 2) * (log(l_aa) + log(l_bb))
  ln_bottom <- (delta_star + numSq2 + 1) / 2 * (log(l_aa ^2 * l_ab ^2 + l_aa ^2 * l_bb ^2))
  ln_gam_bit <- -log(2 * sqrt(pi)) + lgamma((delta_star + numSq2 + 1) / 2 ) - lgamma((delta_star + numSq2 + 2) / 2)
  ln_h <- ln_top - ln_bottom + ln_gam_bit
  output <- ln_h
  
  return(output)
}

########## generate_graph_zo ##########

generate_graph_zo <- function(g_current, delta, Phi, N, S_N, warmup, iters_to_keep)
{
  p <- max(dim(Phi)[1], dim(Phi)[2])
  iter <- iters_to_keep
  if (iter <= warmup)
  {
    stop('not enough iterations specified')
  }
  sample <- iter - warmup
  counter_for_keeping <- 0
  graph_cumulative <- matrix(0, p, p)
  graph_pm <- matrix(0, p, p)
  graph_iter <- array(0, c(p, p, iters_to_keep))
  graph_tally <- array(0, c(p, p, 1))
  num_of_diff_graphs <- 1
  accept_count_warmup <- 0
  accept_count_sample <- 0
  accept_percent_warmup <- 0
  accept_percent_sample <- 0
  g_next <- g_current
  for (i in 1 : iter)
  {
    equal_indicator <- 0
    g_current <- g_next
    ss <- check_chordal(g_current)
    chordal <- ss[[1]]
    order <- ss[[2]]
    cliques <- chordal_to_ripcliques_zo(g_current, order)
    jtree <- ripcliques_to_jtree_zo(cliques)
    sepsize <- sepsize_zo(cliques, jtree)
    reach_graph <- reachability_graph(g_current)
    u <- runif(1)
    accept_prob <- 0
    output <- next_graph_candidate_zo(g_current, jtree, sepsize, cliques, reach_graph)
    g_proposal <- output[[1]]
    a <- output[[2]]
    b <- output[[3]]
    C <- output[[4]]
    CASE_add <- output[[5]]
    CASE_delete <- output[[6]]
    if (CASE_add == 1)
    {
      ln_likelihood_ratio <- ln_h_ratio_zo(C, a, b, delta, Phi) - ln_h_ratio_zo(C, a, b, delta + N, Phi + S_N)
      likelihood_ratio <- exp(ln_likelihood_ratio)
    }else{
      if (CASE_delete == 1)
      {
        ln_likelihood_ratio <- ln_h_ratio_zo(C, a, b, delta + N, Phi + S_N) - ln_h_ratio_zo(C, a, b, delta, Phi)
        likelihood_ratio <- exp(ln_likelihood_ratio)
      }
    }
    accept_prob <- min(1, likelihood_ratio)
    if (u <= accept_prob)
    {
      g_next <- g_proposal
      if (i <= warmup)
      {
        accept_count_warmup <- accept_count_warmup + 1
      }else{
        accept_count_sample <- accept_count_sample + 1
      }
    }else{
      g_next <- g_current
    }
  }
  
  accept_percent_sample <- accept_count_sample / iters_to_keep
  output <- list(g_next, accept_percent_sample)
  
  return(output)
}


########## Sim_HIW ##########

Sim_HIW <- function(G, S, n)
{
  TT <- max(dim(G)[1], dim(G)[2])
  cc <- check_chordal(G)
  check <- cc[[1]]
  order <- cc[[2]]
  if (check == 0)
  {
    stop('graph not decomp')
  }
  cliques <- chordal_to_ripcliques_zo(G, order)
  Tost_Thi_hat <- g_constrain_zo(S, order, cliques)[[1]]
  sigma_identity <- generate_HIW_g_delta_identity_zo(G, cliques, n)[[1]]
  mysigma <- transform_g_conditional_HIW_no_mcs_zo(sigma_identity, G, cliques, n, diag(rep(1, TT)), Tost_Thi_hat)
  output <- list(mysigma[[1]], mysigma[[2]])
  
  return(output)
}

########## setdiff_zo ##########

setdiff_zo <- function(a1, a2)
{
  b <- a1 - a2
  if (length(b) > 0)
  {
    for (j in 1 : length(b))
    {
      if (b[j] < 0)
      {
        b[j] <- a1[j]
      }
    }
  }
  output <- b
  
  return(output)
}

########## setdiag ##########

setdiag <- function(M, v)
{
  diag(M) <- v
  output <- M
  
  return(output)
}

########## neighbours_vertex_cell ##########

neighbours_vertex_cell <- function(G, i)
{
  
  return(which(G[i, ] != 0))
}

########## union_zo ##########

union_zo <- function(v1, v2)
{
  b <- v1 + v2
  if (length(b) > 0)
  {
    for (j in 1 : length(b))
    {
      if (b[j] > 1) 
      {
        b[j] <- 1
      }
    }
  }
  output <- b
  
  return(output)
}

########## neighbours_vertex_zo ##########

neighbours_vertex_zo <- function(G, vi)
{
  
  return(G[, vi])
}

########## intersect_zo ##########

intersect_zo <- function(a1, a2)
{
  
  return(a1 * a2)
}

########## check_chordal ##########

check_chordal <- function(g)
{
  g <- setdiag(g, 1)
  p <- dim(g)[1]
  order <- rep(0, p)
  chordal <- 1
  numbered <- 1
  order[1] <- 1
  for (i in 2 : p)
  {
    capU <- setdiff(1 : p, numbered)
    score <- rep(0, length(capU))
    if (length(capU) > 0)
    {
      for (u_i in 1 : length(capU))
      {
        u <- capU[u_i]
        score[u_i] <- length(intersect(neighbours_vertex_cell(g, u), numbered))
      }
    }
    argmax <- which(score == max(score))[1]
    u <- capU[argmax]
    numbered <- c(numbered, u)
    order[i] <- u
    pa <- intersect(neighbours_vertex_cell(g, u), order[1 : (i - 1)])
    if (!identical(matrix(g[pa, pa], length(pa), length(pa)), matrix(1, length(pa), length(pa))))
    {
      chordal <- 0
      
      break
    }
  }
  output <- list(chordal, order)
  
  return(output)
}

########## chordal_to_ripcliques_zo ##########

chordal_to_ripcliques_zo <- function(G, order)
{
  p <- dim(G)[1]
  pa <- matrix(0, p, p)
  num_pa <- rep(0, p)
  for (i in 2 : p)
  {
    v <- order[i]
    pre_v <- rep(0, p)
    pre_v[order[1 : (i - 1)]] <- 1
    ns <- neighbours_vertex_zo(G, v)
    pa_i <- intersect_zo(ns, pre_v)
    num_pa[i] <- length(which(pa_i != 0))
    pa[, v] <- pa_i
  }
  ladder <- rep(0, p)
  for (i in 1 : p)
  {
    if ((i == p) | (sum(pa[, order[i]]) >= sum(pa[, order[i + 1]])))
    {
      ladder[i] <- order[i]
    }
  }
  cliques <- matrix(0, p, p)
  index_non_zero_ladder <- 0
  for (i in 1 : p)
  {
    if (ladder[i] != 0)
    {
      v_col_i <- rep(0, p)
      v_col_i[ladder[i]] <- 1
      index_non_zero_ladder <- index_non_zero_ladder + 1
      cliques[, index_non_zero_ladder] <- union_zo(v_col_i, pa[, ladder[i]])
    }
  }
  output <- cliques
  
  return(output)
}

########## seps_resids_hists_zo ##########

seps_resids_hists_zo <- function(cliques)
{
  p <- dim(cliques)[1]
  num_cliques <- max(which(colSums(cliques) !=0))
  seps <- matrix(0, p, p)
  resids <- hists <- matrix(0, p, p)
  hists[, 1] <- cliques[, 1]
  if (num_cliques >= 2)
  {
    for (index in 2 : num_cliques)
    {
      hists[, index] <- union_zo((cliques[, index]), (hists[, index - 1]))
      seps[, index] <- intersect_zo((cliques[, index]), (hists[, index - 1]))
      resids[, index] <- setdiff_zo((cliques[, index]), (hists[, index - 1]))
    }
  }
  output <- list(seps, resids, hists)
  
  return(output)
}

########## generate_HIW_g_delta_identity_zo ##########

generate_HIW_g_delta_identity_zo <- function(g, cliques, delta)
{
  index_finish <- 0
  index_start <- 0
  num_cliques <- max(which(colSums(cliques) != 0))
  ss <- seps_resids_hists_zo(cliques)
  seps <- ss[[1]]
  residuals <- ss[[2]]
  histories <- ss[[3]]
  num_Rj <- rep(0, num_cliques)
  p <- max(dim(g)[1], dim(g)[2])
  perfect_order <- rep(0, p)
  c_1 <- which(cliques[, 1] != 0)
  perfect_order[1 : length(c_1)] <- c_1
  index_finish <- length(c_1)
  if (num_cliques >= 2)
  {
    for (j in 2 : num_cliques)
    {
      Rj <- which(residuals[, j] !=0 )
      num_Rj[j] <- length(Rj)
      index_start <- index_finish + 1
      index_finish <- index_start + num_Rj[j] - 1
      perfect_order[index_start : index_finish] <- Rj
    }
  }
  rev_perf <- rep(0, p)
  index <- 0
  for (i in 1 : p)
  {
    rev_perf[i] <- perfect_order[p - index]
    index <- index + 1
  }
  g_rev_perf <- matrix(0, p, p)
  g_rev_perf <- g[rev_perf, rev_perf]
  Psi <- matrix(0, p, p)
  for (i in 1 : p)
  {
    if (i == p)
    {
      Psi[p, p] <- rchisq(1, df = delta) ^0.5
    }else{
      nu_i <- sum(g_rev_perf[i, (i + 1) : p])
      Psi[i, i] <- rchisq(1, df = delta + nu_i) ^0.5
    }
  }
  for (i in 1 : (p - 1))
  {
    for (j in (i + 1) : p)
    {
      if (g_rev_perf[i, j] == 0)
      {
        Psi[i, j] <- 0
      }else{
        Psi[i, j] <- rnorm(1)
      }
    }
  }
  K_rev_perf <- sigma_id_rev_perf <- matrix(0, p, p)
  K_rev_perf <- t(Psi) %*% Psi
  sigma_id_rev_perf <- solve(K_rev_perf)
  inverse_permute <- rep(0, p)
  for (j in 1 : p)
  {
    inverse_permute[j] <- which(rev_perf == j)
  }
  K_id <- sigma_id <- matrix(0, p, p)
  K_id <- K_rev_perf[inverse_permute, inverse_permute]
  sigma_id <- sigma_id_rev_perf[inverse_permute, inverse_permute]
  output <- list(sigma_id, K_id, sigma_id_rev_perf, K_rev_perf)
  
  return(output)
}

########## g_constrain_zo ##########

g_constrain_zo <- function(B, order, cliques)
{
  p <- dim(B)[1]
  seps <- seps_resids_hists_zo(cliques)[[1]]
  num_cliques <- max(which(colSums(cliques) != 0))
  sum_big_K_cliques <- sum_big_K_seps <- matrix(0, p, p)
  if (num_cliques >= 1)
  {
    for (j in 1 : num_cliques)
    {
      cj <- which(cliques[, j] != 0)
      B_cj <- B[cj, cj]
      K_cj <- solve(B_cj)
      big_Kj <- matrix(0, p, p)
      big_Kj[cj, cj] <- K_cj
      sum_big_K_cliques <- sum_big_K_cliques + big_Kj
    }
  }
  if (num_cliques >= 1)
  {
    for (j in 1 : num_cliques)
    {
      sj <- which(seps[, j] != 0)
      B_sj <- B[sj, sj]
      if (length(sj) > 0)
      {
        K_sj <- solve(B_sj)
      }
      big_Kj <- matrix(0, p, p)
      if (length(sj) > 0)
      {
        big_Kj[sj, sj] <- K_sj
      }
      sum_big_K_seps <- sum_big_K_seps + big_Kj
    }
  }
  K_hat <- sum_big_K_cliques - sum_big_K_seps
  B_hat <- solve(K_hat + diag(rep(0.00000001, dim(K_hat)[2])))
  output <- list(B_hat, K_hat)
  
  return(output)
}

########## transform_g_conditional_HIW_no_mcs_zo ##########

transform_g_conditional_HIW_no_mcs_zo <- function(sigma_B, g, cliques, delta, B, D)
{
  p <- dim(g)[1]
  num_cliques <- max(which(colSums(cliques) != 0))
  ss <- seps_resids_hists_zo(cliques)
  seps <- ss[[1]]
  residuals <- ss[[2]]
  histories <- ss[[3]]
  num_Cj <- rep(0, num_cliques)
  num_Sj <- rep(0, num_cliques)
  num_Rj <- rep(0, num_cliques)
  num_Cj[1] <- sum(cliques[, 1])
  perfect_order <- rep(0, p)
  perfect_order[1 : num_Cj[1]] <- which(cliques[, 1] !=0 )
  index_finish <- num_Cj[1]
  if (num_cliques >= 2)
  {
    for (j in 2 : num_cliques)
    {
      Cj <- which(cliques[, j] != 0)
      Rj <- which(residuals[, j] != 0)
      Sj <- which(seps[, j] != 0)
      num_Cj[j] <- length(Cj)
      num_Rj[j] <- length(Rj)
      num_Sj[j] <- length(Sj)
      index_start <- index_finish + 1
      index_finish <- index_start + num_Rj[j] - 1
      perfect_order[index_start : index_finish] <- Rj
    }
  }
  rev_perf <- rep(0, p)
  index <- 0
  for (i in 1 : p)
  {
    rev_perf[i] <- perfect_order[p - index]
    index <- index + 1
  }
  g_rev_perf <- g[rev_perf, rev_perf]
  B_rev_perf <- B[rev_perf, rev_perf]
  D_rev_perf <- D[rev_perf, rev_perf]
  sigma_B_rev_perf <- sigma_B[rev_perf, rev_perf]
  K_B <- solve(sigma_B_rev_perf)
  choleskyK_B <- chol(K_B)
  c1 <- num_Cj[1]
  index_start <- p - c1 + 1
  index_finish <- p
  indexC1_in_Upsilon_D <- rep(0, c1)
  indexC1_in_Upsilon_D <- index_start : index_finish
  Upsilon_D <- matrix(0, p, p)
  B_1 <- D_1 <- Q_1 <- P_1 <- O_1 <- matrix(0, c1, c1)
  B_1 <- B_rev_perf[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D]
  Q_1 <- chol(solve(B_1))
  D_1 <- D_rev_perf[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D]
  P_1 <- chol(solve(D_1))
  O_1 <- solve(Q_1) %*% P_1
  Upsilon_D[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D] <- choleskyK_B[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D] %*% O_1
  if (num_cliques >= 2)
  {
    for (j in 2 : num_cliques)
    {
      cj <- num_Cj[j]
      indexCj_in_Upsilon_D <- rep(0, cj)
      Rj <- which(residuals[, j] != 0)
      rj <- num_Rj[j]
      indexRj_in_Upsilon_D <- rep(0, rj)
      indexRj_inOj <- rep(0, rj)
      unsort_indexRj_in_Upsilon_D <- rep(0, rj)
      Sj <- which(seps[, j] != 0)
      sj <- num_Sj[j]
      unsort_indexSj_in_Upsilon_D <- rep(0, sj)
      indexSj_in_Upsilon_D <- rep(0, sj)
      B_j <- D_j <- Q_j <- P_j <- O_j <- matrix(0, cj, cj)
      if (rj > 0)
      {
        for (k in 1 : rj)
        {
          unsort_indexRj_in_Upsilon_D[k] <- which(rev_perf == Rj[k])
        }
        indexRj_in_Upsilon_D <- sort(unsort_indexRj_in_Upsilon_D)
      }
      if (sj > 0)
      {
        for (k in 1 : sj)
        {
          unsort_indexSj_in_Upsilon_D[k] <- which(rev_perf == Sj[k])
        }
      indexSj_in_Upsilon_D <- sort(unsort_indexSj_in_Upsilon_D)
      }
    indexCj_in_Upsilon_D <- c(indexRj_in_Upsilon_D, indexSj_in_Upsilon_D)
    B_j <- (B_rev_perf[indexCj_in_Upsilon_D, indexCj_in_Upsilon_D])
    Q_j <- chol(solve(B_j))
    D_j <- D_rev_perf[indexCj_in_Upsilon_D, indexCj_in_Upsilon_D]
    P_j <- chol(solve(D_j))
    O_j <- solve(Q_j) %*% P_j
    if (rj > 0)
    {
      indexRj_inOj <- 1 : rj
    }else{
      indexRj_inOj <- NULL
    }
    if ((rj + 1) <= cj)
    {
      indexSj_inOj <- (rj + 1) : cj
    }else{
      indexSj_inOj <- NULL
    }
    Upsilon_D[indexRj_in_Upsilon_D, indexRj_in_Upsilon_D] <- as.matrix(choleskyK_B[indexRj_in_Upsilon_D, indexRj_in_Upsilon_D]) %*% as.matrix(O_j[indexRj_inOj, indexRj_inOj])
    Upsilon_D[indexRj_in_Upsilon_D, indexSj_in_Upsilon_D] <- as.matrix(choleskyK_B[indexRj_in_Upsilon_D, indexRj_in_Upsilon_D]) %*% O_j[indexRj_inOj,indexSj_inOj] + choleskyK_B[indexRj_in_Upsilon_D, indexSj_in_Upsilon_D] %*% as.matrix(O_j[indexSj_inOj, indexSj_inOj])
    }
  }
  K_D_rev_perf <- sigma_D_rev_perf <- matrix(0, p, p)
  K_D_rev_perf <- t(Upsilon_D) %*% Upsilon_D
  sigma_D_rev_perf <- solve(K_D_rev_perf)
  inverse_permute <- rep(0, p)
  for (j in 1 : p)
  {
    inverse_permute[j] <- which(rev_perf == j)
  }
  K_D <- sigma_D <- matrix(0, p, p)
  K_D <- K_D_rev_perf[inverse_permute, inverse_permute]
  sigma_D <- sigma_D_rev_perf[inverse_permute, inverse_permute]
  output <- list(sigma_D, K_D)
  
  return(output)
}


########## Sim_HIW ##########

Sim_HIW <- function(G, S, n)
{
  TT <- max(dim(G)[1], dim(G)[2])
  cc <- check_chordal(G)
  check <- cc[[1]]
  order <- cc[[2]]
  if (check == 0)
  {
    stop('graph not decomp')
  }
  cliques <- chordal_to_ripcliques_zo(G, order)
  Tost_Thi_hat <- g_constrain_zo(S, order, cliques)[[1]]
  sigma_identity <- generate_HIW_g_delta_identity_zo(G, cliques, n)[[1]]
  mysigma <- transform_g_conditional_HIW_no_mcs_zo(sigma_identity, G, cliques, n, diag(rep(1, TT)), Tost_Thi_hat)
  output <- list(mysigma[[1]], mysigma[[2]])
  
  return(output)
}

########## setdiff_zo ##########

setdiff_zo <- function(a1, a2)
{
  b <- a1 - a2
  if (length(b) > 0)
  {
    for (j in 1 : length(b))
    {
      if (b[j] < 0)
      {
        b[j] <- a1[j]
      }
    }
  }
  output <- b
  
  return(output)
}

########## setdiag ##########

setdiag <- function(M, v)
{
  diag(M) <- v
  output <- M
  
  return(output)
}

########## neighbours_vertex_cell ##########

neighbours_vertex_cell <- function(G, i)
{
  
  return(which(G[i, ] != 0))
}

########## union_zo ##########

union_zo <- function(v1, v2)
{
  b <- v1 + v2
  if (length(b) > 0)
  {
    for (j in 1 : length(b))
    {
      if (b[j] > 1) 
      {
        b[j] <- 1
      }
    }
  }
  output <- b
  
  return(output)
}

########## neighbours_vertex_zo ##########

neighbours_vertex_zo <- function(G, vi)
{
  
  return(G[, vi])
}

########## intersect_zo ##########

intersect_zo <- function(a1, a2)
{
  
  return(a1 * a2)
}

########## check_chordal ##########

check_chordal <- function(g)
{
  g <- setdiag(g, 1)
  p <- dim(g)[1]
  order <- rep(0, p)
  chordal <- 1
  numbered <- 1
  order[1] <- 1
  for (i in 2 : p)
  {
    capU <- setdiff(1 : p, numbered)
    score <- rep(0, length(capU))
    if (length(capU) > 0)
    {
      for (u_i in 1 : length(capU))
      {
        u <- capU[u_i]
        score[u_i] <- length(intersect(neighbours_vertex_cell(g, u), numbered))
      }
    }
    argmax <- which(score == max(score))[1]
    u <- capU[argmax]
    numbered <- c(numbered, u)
    order[i] <- u
    pa <- intersect(neighbours_vertex_cell(g, u), order[1 : (i - 1)])
    if (!identical(matrix(g[pa, pa], length(pa), length(pa)), matrix(1, length(pa), length(pa))))
    {
      chordal <- 0
      
      break
    }
  }
  output <- list(chordal, order)
  
  return(output)
}

########## chordal_to_ripcliques_zo ##########

chordal_to_ripcliques_zo <- function(G, order)
{
  p <- dim(G)[1]
  pa <- matrix(0, p, p)
  num_pa <- rep(0, p)
  for (i in 2 : p)
  {
    v <- order[i]
    pre_v <- rep(0, p)
    pre_v[order[1 : (i - 1)]] <- 1
    ns <- neighbours_vertex_zo(G, v)
    pa_i <- intersect_zo(ns, pre_v)
    num_pa[i] <- length(which(pa_i != 0))
    pa[, v] <- pa_i
  }
  ladder <- rep(0, p)
  for (i in 1 : p)
  {
    if ((i == p) | (sum(pa[, order[i]]) >= sum(pa[, order[i + 1]])))
    {
      ladder[i] <- order[i]
    }
  }
  cliques <- matrix(0, p, p)
  index_non_zero_ladder <- 0
  for (i in 1 : p)
  {
    if (ladder[i] != 0)
    {
      v_col_i <- rep(0, p)
      v_col_i[ladder[i]] <- 1
      index_non_zero_ladder <- index_non_zero_ladder + 1
      cliques[, index_non_zero_ladder] <- union_zo(v_col_i, pa[, ladder[i]])
    }
  }
  output <- cliques
  
  return(output)
}

########## seps_resids_hists_zo ##########

seps_resids_hists_zo <- function(cliques)
{
  p <- dim(cliques)[1]
  num_cliques <- max(which(colSums(cliques) !=0))
  seps <- matrix(0, p, p)
  resids <- hists <- matrix(0, p, p)
  hists[, 1] <- cliques[, 1]
  if (num_cliques >= 2)
  {
    for (index in 2 : num_cliques)
    {
      hists[, index] <- union_zo((cliques[, index]), (hists[, index - 1]))
      seps[, index] <- intersect_zo((cliques[, index]), (hists[, index - 1]))
      resids[, index] <- setdiff_zo((cliques[, index]), (hists[, index - 1]))
    }
  }
  output <- list(seps, resids, hists)
  
  return(output)
}

########## generate_HIW_g_delta_identity_zo ##########

generate_HIW_g_delta_identity_zo <- function(g, cliques, delta)
{
  index_finish <- 0
  index_start <- 0
  num_cliques <- max(which(colSums(cliques) != 0))
  ss <- seps_resids_hists_zo(cliques)
  seps <- ss[[1]]
  residuals <- ss[[2]]
  histories <- ss[[3]]
  num_Rj <- rep(0, num_cliques)
  p <- max(dim(g)[1], dim(g)[2])
  perfect_order <- rep(0, p)
  c_1 <- which(cliques[, 1] != 0)
  perfect_order[1 : length(c_1)] <- c_1
  index_finish <- length(c_1)
  if (num_cliques >= 2)
  {
    for (j in 2 : num_cliques)
    {
      Rj <- which(residuals[, j] !=0 )
      num_Rj[j] <- length(Rj)
      index_start <- index_finish + 1
      index_finish <- index_start + num_Rj[j] - 1
      perfect_order[index_start : index_finish] <- Rj
    }
  }
  rev_perf <- rep(0, p)
  index <- 0
  for (i in 1 : p)
  {
    rev_perf[i] <- perfect_order[p - index]
    index <- index + 1
  }
  g_rev_perf <- matrix(0, p, p)
  g_rev_perf <- g[rev_perf, rev_perf]
  Psi <- matrix(0, p, p)
  for (i in 1 : p)
  {
    if (i == p)
    {
      Psi[p, p] <- rchisq(1, df = delta) ^0.5
    }else{
      nu_i <- sum(g_rev_perf[i, (i + 1) : p])
      Psi[i, i] <- rchisq(1, df = delta + nu_i) ^0.5
    }
  }
  for (i in 1 : (p - 1))
  {
    for (j in (i + 1) : p)
    {
      if (g_rev_perf[i, j] == 0)
      {
        Psi[i, j] <- 0
      }else{
        Psi[i, j] <- rnorm(1)
      }
    }
  }
  K_rev_perf <- sigma_id_rev_perf <- matrix(0, p, p)
  K_rev_perf <- t(Psi) %*% Psi
  sigma_id_rev_perf <- solve(K_rev_perf)
  inverse_permute <- rep(0, p)
  for (j in 1 : p)
  {
    inverse_permute[j] <- which(rev_perf == j)
  }
  K_id <- sigma_id <- matrix(0, p, p)
  K_id <- K_rev_perf[inverse_permute, inverse_permute]
  sigma_id <- sigma_id_rev_perf[inverse_permute, inverse_permute]
  output <- list(sigma_id, K_id, sigma_id_rev_perf, K_rev_perf)
  
  return(output)
}

########## g_constrain_zo ##########

g_constrain_zo <- function(B, order, cliques)
{
  p <- dim(B)[1]
  seps <- seps_resids_hists_zo(cliques)[[1]]
  num_cliques <- max(which(colSums(cliques) != 0))
  sum_big_K_cliques <- sum_big_K_seps <- matrix(0, p, p)
  if (num_cliques >= 1)
  {
    for (j in 1 : num_cliques)
    {
      cj <- which(cliques[, j] != 0)
      B_cj <- B[cj, cj]
      K_cj <- solve(B_cj)
      big_Kj <- matrix(0, p, p)
      big_Kj[cj, cj] <- K_cj
      sum_big_K_cliques <- sum_big_K_cliques + big_Kj
    }
  }
  if (num_cliques >= 1)
  {
    for (j in 1 : num_cliques)
    {
      sj <- which(seps[, j] != 0)
      B_sj <- B[sj, sj]
      if (length(sj) > 0)
      {
        K_sj <- solve(B_sj)
      }
      big_Kj <- matrix(0, p, p)
      if (length(sj) > 0)
      {
        big_Kj[sj, sj] <- K_sj
      }
      sum_big_K_seps <- sum_big_K_seps + big_Kj
    }
  }
  K_hat <- sum_big_K_cliques - sum_big_K_seps
  B_hat <- solve(K_hat + diag(rep(0.00000001, dim(K_hat)[2])))
  output <- list(B_hat, K_hat)
  
  return(output)
}

########## transform_g_conditional_HIW_no_mcs_zo ##########

transform_g_conditional_HIW_no_mcs_zo <- function(sigma_B, g, cliques, delta, B, D)
{
  p <- dim(g)[1]
  num_cliques <- max(which(colSums(cliques) != 0))
  ss <- seps_resids_hists_zo(cliques)
  seps <- ss[[1]]
  residuals <- ss[[2]]
  histories <- ss[[3]]
  num_Cj <- rep(0, num_cliques)
  num_Sj <- rep(0, num_cliques)
  num_Rj <- rep(0, num_cliques)
  num_Cj[1] <- sum(cliques[, 1])
  perfect_order <- rep(0, p)
  perfect_order[1 : num_Cj[1]] <- which(cliques[, 1] !=0 )
  index_finish <- num_Cj[1]
  if (num_cliques >= 2)
  {
    for (j in 2 : num_cliques)
    {
      Cj <- which(cliques[, j] != 0)
      Rj <- which(residuals[, j] != 0)
      Sj <- which(seps[, j] != 0)
      num_Cj[j] <- length(Cj)
      num_Rj[j] <- length(Rj)
      num_Sj[j] <- length(Sj)
      index_start <- index_finish + 1
      index_finish <- index_start + num_Rj[j] - 1
      perfect_order[index_start : index_finish] <- Rj
    }
  }
  rev_perf <- rep(0, p)
  index <- 0
  for (i in 1 : p)
  {
    rev_perf[i] <- perfect_order[p - index]
    index <- index + 1
  }
  g_rev_perf <- g[rev_perf, rev_perf]
  B_rev_perf <- B[rev_perf, rev_perf]
  D_rev_perf <- D[rev_perf, rev_perf]
  sigma_B_rev_perf <- sigma_B[rev_perf, rev_perf]
  K_B <- solve(sigma_B_rev_perf)
  choleskyK_B <- chol(K_B)
  c1 <- num_Cj[1]
  index_start <- p - c1 + 1
  index_finish <- p
  indexC1_in_Upsilon_D <- rep(0, c1)
  indexC1_in_Upsilon_D <- index_start : index_finish
  Upsilon_D <- matrix(0, p, p)
  B_1 <- D_1 <- Q_1 <- P_1 <- O_1 <- matrix(0, c1, c1)
  B_1 <- B_rev_perf[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D]
  Q_1 <- chol(solve(B_1))
  D_1 <- D_rev_perf[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D]
  P_1 <- chol(solve(D_1))
  O_1 <- solve(Q_1) %*% P_1
  Upsilon_D[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D] <- choleskyK_B[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D] %*% O_1
  if (num_cliques >= 2)
  {
    for (j in 2 : num_cliques)
    {
      cj <- num_Cj[j]
      indexCj_in_Upsilon_D <- rep(0, cj)
      Rj <- which(residuals[, j] != 0)
      rj <- num_Rj[j]
      indexRj_in_Upsilon_D <- rep(0, rj)
      indexRj_inOj <- rep(0, rj)
      unsort_indexRj_in_Upsilon_D <- rep(0, rj)
      Sj <- which(seps[, j] != 0)
      sj <- num_Sj[j]
      unsort_indexSj_in_Upsilon_D <- rep(0, sj)
      indexSj_in_Upsilon_D <- rep(0, sj)
      B_j <- D_j <- Q_j <- P_j <- O_j <- matrix(0, cj, cj)
      if (rj > 0)
      {
        for (k in 1 : rj)
        {
          unsort_indexRj_in_Upsilon_D[k] <- which(rev_perf == Rj[k])
        }
        indexRj_in_Upsilon_D <- sort(unsort_indexRj_in_Upsilon_D)
      }
      if (sj > 0)
      {
        for (k in 1 : sj)
        {
          unsort_indexSj_in_Upsilon_D[k] <- which(rev_perf == Sj[k])
        }
      indexSj_in_Upsilon_D <- sort(unsort_indexSj_in_Upsilon_D)
      }
    indexCj_in_Upsilon_D <- c(indexRj_in_Upsilon_D, indexSj_in_Upsilon_D)
    B_j <- (B_rev_perf[indexCj_in_Upsilon_D, indexCj_in_Upsilon_D])
    Q_j <- chol(solve(B_j))
    D_j <- D_rev_perf[indexCj_in_Upsilon_D, indexCj_in_Upsilon_D]
    P_j <- chol(solve(D_j))
    O_j <- solve(Q_j) %*% P_j
    if (rj > 0)
    {
      indexRj_inOj <- 1 : rj
    }else{
      indexRj_inOj <- NULL
    }
    if ((rj + 1) <= cj)
    {
      indexSj_inOj <- (rj + 1) : cj
    }else{
      indexSj_inOj <- NULL
    }
    Upsilon_D[indexRj_in_Upsilon_D, indexRj_in_Upsilon_D] <- as.matrix(choleskyK_B[indexRj_in_Upsilon_D, indexRj_in_Upsilon_D]) %*% as.matrix(O_j[indexRj_inOj, indexRj_inOj])
    Upsilon_D[indexRj_in_Upsilon_D, indexSj_in_Upsilon_D] <- as.matrix(choleskyK_B[indexRj_in_Upsilon_D, indexRj_in_Upsilon_D]) %*% O_j[indexRj_inOj,indexSj_inOj] + choleskyK_B[indexRj_in_Upsilon_D, indexSj_in_Upsilon_D] %*% as.matrix(O_j[indexSj_inOj, indexSj_inOj])
    }
  }
  K_D_rev_perf <- sigma_D_rev_perf <- matrix(0, p, p)
  K_D_rev_perf <- t(Upsilon_D) %*% Upsilon_D
  sigma_D_rev_perf <- solve(K_D_rev_perf)
  inverse_permute <- rep(0, p)
  for (j in 1 : p)
  {
    inverse_permute[j] <- which(rev_perf == j)
  }
  K_D <- sigma_D <- matrix(0, p, p)
  K_D <- K_D_rev_perf[inverse_permute, inverse_permute]
  sigma_D <- sigma_D_rev_perf[inverse_permute, inverse_permute]
  output <- list(sigma_D, K_D)
  
  return(output)
}


########## mod_BDgraphInit ##########

mod_BDgraphInit <- function(data, samplerPar, 
                            n = NULL, not.cont = NULL, df.prior = 2, g.prior = 0.5, 
                            threshold = 1e-08)
{
  burnin <- samplerPar$BDgraphPar$burnin
  niter <- samplerPar$BDgraphPar$niter
  cores <- samplerPar$BDgraphPar$cores
  
  if (niter < burnin)
  {
    stop("Number of iterations must be higher than the burn-in")
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
  G <- BDgraph::get_g_start(g.start = "empty", g_prior = g_prior, p = p)
  G[g_prior == 1] <- 1
  G[g_prior == 0] <- 0
  G[lower.tri(G, diag(TRUE))] <- 0
  G <- G + t(G)
  K <- BDgraph::get_K_start(G = G, g.start = "empty", Ts = Ts, b_star = b_star, threshold = threshold)
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
