# 135 characters #####################################################################################################################
#' @title MR2: Multi-response Mendelian randomization
#' @description Markov chain Monte Carlo implementation of Bayesian Multi-response Mendelian randomization model for two-sample 
#' summary-level data
#'
#' @param Y (number of IVs) times (number of responses) matrix of summary-level responses
#' @param X (number of IVs) times (number of exposures) matrix of summary-level exposures
#' @param EVgamma A priori number of expected exposures directly associated with each response and its variance. For details, see 
#' \insertCite{Kohn2001;textual}{MR2}, \insertCite{Alexopoulos2021;textual}{MR2} and \insertCite{Zuber2023;textual}{MR2}. Default 
#' value is \code{EVgamma = c(2, 2)} which implies a priori range of important exposures associated with each response between 0 and 6
#' @param tau Prior variance of the (direct causal effect) theta prior with the default value set at \code{1}. If \code{std = TRUE}, 
#' then \code{tau = 1}
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
#' @param std Logical parameter to standardise the exposures before the analysis. Default value is set at \code{TRUE}
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
#'   \item{\code{postMean}}{ List of the posterior means of \code{list(theta, G, R, D)} }
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
#' # Example 1: Analysis of one replication of simulated Scenario II-Counfounding with q = 5 
#' # responses, p = 15 exposures and nIV = 100 genetic variants used as IVs. The number of expected 
#' # exposures directly associated with each response is set at 2 and its variance at 2, with a 
#' # priori range of direct causal associations ranging between 0 and 8
#'
#' Sim_Counfounding <- Scenario_Counfounding(nRep = 1, q = 5, p = 15, nIV = 100, seed = 280610971)
#' beta_Y <- Sim_Counfounding$Y
#' beta_X <- Sim_Counfounding$X
#'
#' MR2_output <- MR2(beta_Y, beta_X, EVgamma = c(2, 2), 
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
#' beta_Y <- Sim_Pleiotropy$Y
#' beta_X <- Sim_Pleiotropy$X
#'
#' MR2_output <- MR2(beta_Y, beta_X, EVgamma = c(1, 3),
#'                   niter = 15000, burnin = 5000, thin = 10, monitor = 1000, 
#'                   nonDecomp = TRUE, seed = 280610971)
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
#'             rnorm runif var
#' @importFrom truncnorm rtruncnorm
#' @importFrom utils capture.output
#' @importFrom Rdpack reprompt


MR2 <- function(Y, X, 
                EVgamma = c(2, 2), tau = 1, 
                niter, burnin, thin, monitor, 
                fullCov = FALSE, nonDecomp = FALSE, 
                alpha = 0.05, probAddDel = 0.9, probAdd = 0.5, probMix = 0.25, geomMean = EVgamma[1], 
                graphParNiter = max(16, dim(Y)[2] * dim(Y)[2]), std = TRUE, seed = 31122021)
{
  set.seed(seed)
  
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  
  if (std == TRUE)
  {
    X <- scale(X, center = FALSE, scale = TRUE)
    tau <- 1
  }
  
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
  samplerPar$link <- link
  samplerPar$alpha <- alpha
  samplerPar$probAddDel <- probAddDel
  samplerPar$probAdd <- probAdd
  samplerPar$probMix <- probMix
  samplerPar$geomMean <- geomMean
  
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
  
  Eqg <- EVgamma[1]
  Varqg <- EVgamma[2]
  Epi <- Eqg / pStar
  betaMean <- Epi
  betaVar <- (Varqg - pStar * Epi * (1 - pStar * Epi)) / (pStar * (pStar - 1) * Epi)
  hyperPar$aPi <- (betaMean * betaVar - betaMean) / (betaMean - betaVar)
  hyperPar$bPi <- hyperPar$aPi * (1 - betaMean) / betaMean
  
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
      nOlvs <- max(Y[, k], na.rm = TRUE) + 1
      nCat[k] <- nOlvs
      cutPoint[[k]] <- rep(NA, nOlvs + 1)
      cutPoint[[k]][c(1, 2, nOlvs + 1)] <- c(-Inf, 0, Inf)
      countCat <- 0
      for(kk in 3 : nOlvs)
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
    BD <- mod_BDgraphInit(WCurr, samplerPar)
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
        BD <- mod_BDgraph(data = WCurr, g.start = BD$last_graph, K.start = BD$last_K, samplerPar)        
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
  
  postMean <- NULL
  postMean$gamma <- matrix(colMeans(abs(BSave) > 0), p, m)
  postMean$theta <- matrix(colMeans(BSave), p, m)
  postMean$G <- apply(GSave, c(1, 2), mean)
  postMean$R <- apply(RSave, c(1, 2), mean)
  postMean$D <- colMeans(t(D))
  
  colnames(postMean$gamma) <- colnames_Y
  rownames(postMean$gamma) <- colnames_X
  colnames(postMean$theta) <- colnames_Y
  rownames(postMean$theta) <- colnames_X
  colnames(postMean$G) <- colnames_Y
  rownames(postMean$G) <- colnames_Y
  colnames(postMean$R) <- colnames_Y
  rownames(postMean$R) <- colnames_Y
  names(postMean$D) <- colnames_Y
  
  opt = list(std = std, seed = seed)
  
  output <- list(theta = BSave, G = GSave, R = RSave, D = DSave, 
                 postMean = postMean, 
                 hyperPar = hyperPar,
                 samplerPar = samplerPar, 
                 opt = opt)
  
  return(output)
}
