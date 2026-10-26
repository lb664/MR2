# 135 characters #####################################################################################################################
#
#
#' @title MR2: Multi-response Mendelian randomisation
#' @description Markov chain Monte Carlo implementation of Bayesian Multi-response Mendelian randomisation model for two-sample
#' summary-level data
#'
#' @param betaHat_Y (number of IVs) times (number of responses) matrix of summary-level responses
#' @param betaHat_X (number of IVs) times (number of exposures) matrix of summary-level exposures
#' @param EVgamma A priori number of expected exposures directly associated with each response and its variance
#' @param tau Prior variance of the direct effects theta. Default \code{NULL}: the prior variance is estimated, one value per
#' response, with the hyperprior \eqn{\tau_k \sim \mathrm{Inverse\text{-}Gamma}(\alpha, \eta)}, \eqn{\eta \sim \mathrm{Gamma}(g, h)}
#' (see \code{tauPar}). A single positive number fixes the prior variance at that value for every response
#' @param tauPar Hyperparameters \code{c(alpha, g, h)} of the hyperprior on \code{tau}, used when \code{tau = NULL}.
#' Default \code{c(alpha = 2, g = 0.2, h = 10)}
#' @param niter Number of Markov chain Monte Carlo iterations (including burn-in)
#' @param burnin Number of Markov chain Monte Carlo iterations to be discarded as burn-in
#' @param thin Store the Markov chain Monte Carlo output at every thin-th iteration
#' @param monitor Refresh the console progress bar every monitored-th iteration
#' @param nonDecomp Logical, \code{TRUE} for a non-decomposable graphical model for the inverse correlation matrix
#' @param alpha Level of significance (\code{0.05} default) used to construct the proposal distribution of \code{gamma}
#' @param probAddDel Probability (\code{0.9} default) of adding/deleting one exposure during MCMC exploration
#' @param probAdd Probability (\code{0.5} default) of adding one exposure during MCMC exploration
#' @param probMix Probability (\code{0.25} default) of a geometric proposal for the exposure index being added/deleted
#' @param geomMean Mean of the geometric proposal distribution. Default is \code{EVgamma[1]}
#' @param graphParNiter Number of MCMC internal iterations for the graphical model. Default \code{max(16, m^2)}
#' @param seed Seed used to initialise the Markov chain Monte Carlo algorithm
#'
#' @details
#' \insertCite{Alexopoulos2021;textual}{MR2}
#'
#' @export
#'
#' @return A list \code{list(theta, G, R, D, postMean, timeMR2, hyperPar, samplerPar, opt)}; see package documentation for details
#' of each element
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom grDevices gray
#' @importFrom graphics abline axis hist lines
#' @importFrom MASS mvrnorm
#' @importFrom mvnfast dmvn
#' @importFrom stats cor cov2cor dbeta dbinom dgeom kmeans lm nlm pgeom qgeom quantile rbinom rchisq rgamma rnorm runif sd var
#' @importFrom utils capture.output setTxtProgressBar txtProgressBar
#' @importFrom Rdpack reprompt


MR2 <- function(betaHat_Y, betaHat_X,
                EVgamma = c(2, 2), tau = NULL, tauPar = c(alpha = 2, g = 0.2, h = 10),
                niter, burnin, thin, monitor,
                nonDecomp = FALSE,
                alpha = 0.05, probAddDel = 0.9, probAdd = 0.5, probMix = 0.25, geomMean = EVgamma[1],
                graphParNiter = NULL,
                seed = 31122021)
{
  cat("\n##\n## Multi-response Mendelian randomisation (MR2)\n## Version: 0.2.0\n##\n")
  cat("## Copyright (C) 2025 L. Bottolo and V. Zuber\n##\n## The model selected is:\n##\n")

  if (nonDecomp) cat("## MR2 with non-decomposable graph for unmeasured pleiotropy\n")
  else cat("## MR2 with decomposable graph for unmeasured pleiotropy\n")
  cat("##\n\n")

  set.seed(seed)
  startTime <- Sys.time()

  Y <- as.matrix(betaHat_Y)
  X <- as.matrix(betaHat_X)
  rm(betaHat_Y, betaHat_X)

  if (is.null(graphParNiter)) graphParNiter <- max(16L, ncol(Y)^2L)

  n <- nrow(Y); m <- ncol(Y); p <- ncol(X)

  if (is.null(colnames(Y))) colnames(Y) <- paste0("R",  seq_len(m))
  if (is.null(rownames(Y))) rownames(Y) <- paste0("IV", seq_len(n))
  if (is.null(colnames(X))) colnames(X) <- paste0("E",  seq_len(p))
  if (is.null(rownames(X))) rownames(X) <- rownames(Y)

  colnames_Y <- colnames(Y); rownames_Y <- rownames(Y)
  colnames_X <- colnames(X)

  hyperPar <- NULL
  if (is.null(tau)) {
    if (length(tauPar) != 3 || any(tauPar <= 0)) stop("`tauPar` must be c(alpha, g, h) with positive entries.")
    hyperPar$tauRandom <- TRUE
    hyperPar$tauAlpha <- tauPar[[1]]; hyperPar$tauG <- tauPar[[2]]; hyperPar$tauH <- tauPar[[3]]
    hyperPar$tau <- rep(1, m)
  } else {
    if (length(tau) != 1 || !is.numeric(tau) || tau <= 0) stop("`tau` must be NULL or a single positive number.")
    hyperPar$tauRandom <- FALSE
    hyperPar$tau <- rep(tau, m)
  }

  samplerPar <- NULL
  samplerPar$niter <- niter
  samplerPar$burnin <- burnin
  samplerPar$thin <- thin
  samplerPar$monitor <- monitor
  samplerPar$nonDecomp <- nonDecomp
  samplerPar$alpha <- alpha
  samplerPar$probAddDel <- probAddDel
  samplerPar$probAdd <- probAdd
  samplerPar$probMix <- probMix

  if (length(EVgamma) == 2) {
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
                          dimnames = list(paste0("iter", seq_len(niter)),
                                          apply(expand.grid(x = colnames_X, y = colnames_Y), 1, paste, collapse = ",")))
  saved_iters <- seq(thin, niter, by = thin)
  saved_iters <- saved_iters[saved_iters > burnin]
  n_save <- length(saved_iters)

  BSave <- matrix(NA, n_save, p * m,
                  dimnames = list(paste0("iter", saved_iters),
                                  apply(expand.grid(x = colnames_X, y = colnames_Y), 1, paste, collapse = ",")))
  GSave <- RSave <- array(NA, c(m, m, n_save),
                          dimnames = list(colnames_Y, colnames_Y,
                                          paste0("iter", saved_iters)))
  DSave <- matrix(NA, m, n_save,
                  dimnames = list(colnames_Y,
                                  paste0("iter", saved_iters)))
  tauSave <- matrix(NA, m, n_save,
                    dimnames = list(colnames_Y,
                                    paste0("iter", saved_iters)))
  etaSave <- rep(NA_real_, n_save)

  ########## Missing values ##########

  missC <- any(is.na(Y))

  NA_Y <- vector("list", m)
  for (k in seq_len(m)) NA_Y[[k]] <- which(!is.na(Y[, k]))

  ########## Computing p-values needed for gamma proposal (Gaussian responses only) ##########

  pStar <- p
  pval  <- matrix(NA, pStar, m)

  for (j in seq_len(pStar))
    for (k in seq_len(m)) {
      fit <- lm(Y[, k] ~ 0 + X[, j])
      pval[j, k] <- summary(fit)$coefficients[1, 4]
    }

  samplerPar$pval <- pval
  samplerPar$mlogpval <- -log10(samplerPar$pval)

  ########## Beta-binomial prior gamma ##########

  if (length(EVgamma) == 2)
  {
    Eqg <- EVgamma[1]; Varqg <- EVgamma[2]
    Epi <- Eqg / pStar
    betaMean <- Epi
    betaVar <- (Varqg - pStar * Epi * (1 - pStar * Epi)) / (pStar * (pStar - 1) * Epi)
    hyperPar$aPi <- (betaMean * betaVar - betaMean) / (betaMean - betaVar)
    hyperPar$bPi <- hyperPar$aPi * (1 - betaMean) / betaMean

    if (hyperPar$aPi < 0) { hyperPar$aPi <- Epi; hyperPar$bPi <- NA }
  } else {
    hyperPar$aPi <- EVgamma
    hyperPar$bPi <- NA
  }

  ########## Initial values gamma ##########

  gammaCurr <- rep(0, p * m)
  for (kk in 0:(m - 1))
  {
    thres <- -log10(samplerPar$alpha / pStar)
    gammaCurr[which(samplerPar$mlogpval[, kk + 1] > thres) + kk * p] <- 1
  }

  ########## Initial values decomposable graph, correlation matrix and Beta ##########

  GCurr <- diag(rep(1, m)); diag(GCurr) <- 0
  RCurr <- diag(rep(1, m))
  invRCurr <- solve(RCurr)
  BCurr <- rep(0, p * m)
  BCurr[which(gammaCurr != 0)] <- rnorm(length(which(gammaCurr != 0)), 0, 0.1)
  BMat <- matrix(BCurr, m, p, byrow = TRUE)
  d <- rep(1, m)
  D <- diag(d)
  XB <- X %*% t(BMat)

  invSigma <- diag(1 / d) %*% invRCurr %*% diag(1 / d)
  Sigma <- solve(invSigma)

  ########## Initial values Z ##########

  ZCurr <- (Y - XB) %*% diag(1 / d)
  if (missC) ZCurr[is.na(ZCurr)] <- rnorm(sum(is.na(ZCurr)))

  WCurr <- ZCurr %*% D

  ########## Initial values non-decomposable graph ##########

  BD <- NULL
  if (nonDecomp)
  {
    capture.output(BD <- mod_BDgraphInit(WCurr, samplerPar))
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
  pb <- txtProgressBar(min = 0, max = niter, style = 3)

  for (iter in seq_len(niter))
  {
    Latent <- BVS(Y, NA_Y, X, XB, gammaCurr, BCurr, ZCurr, RCurr, invRCurr, d, hyperPar, samplerPar)

    gammaProp <- Latent[[1]]
    gammaCurr <- Latent[[2]]
    BCurr     <- Latent[[3]]
    ZCurr     <- Latent[[4]]
    XB        <- Latent[[5]]
    hyperPar$tau <- Latent[[6]]
    eta       <- Latent[[7]]

    WCurr <- ZCurr %*% D

    S <- crossprod(WCurr)

    if (nonDecomp) {
      capture.output(
      BD <- mod_BDgraph(data = WCurr, g.start = BD$last_graph, K.start = BD$last_K, samplerPar))
      invSigma <- BD$last_K
      Sigma <- solve(invSigma + diag(rep(10^(-10)), m))
      GCurr <- BD$last_graph
    } else {
      GCurr <- .Sample_G(WCurr, GCurr, samplerPar)
      HIW <- Sim_HIW(GCurr, S + diag(rep(1, m)), n + 2)
      Sigma <- HIW[[1]]
      invSigma <- HIW[[2]]
    }

    d <- sqrt(diag(Sigma))
    RCurr <- diag(1 / d) %*% Sigma %*% diag(1 / d)
    invRCurr <- diag(d) %*% invSigma %*% diag(d)
    D <- diag(d)
    ZCurr <- WCurr %*% diag(1 / d)

    gammaPropSave[iter, ] <- gammaProp

    if (iter > burnin && iter %% thin == 0)
    {
      countIter <- countIter + 1
      BSave[countIter, ] <- BCurr
      GSave[, , countIter] <- GCurr
      RSave[, , countIter] <- RCurr
      DSave[, countIter] <- d
      tauSave[, countIter] <- hyperPar$tau
      etaSave[countIter] <- eta
    }
    if (iter %% monitor == 0) setTxtProgressBar(pb, iter)
  }

  setTxtProgressBar(pb, niter)
  close(pb)
  cat("\n")

  DSave <- t(DSave)
  tauSave <- t(tauSave)

  ########## Gamma proposal parameters ##########

  samplerPar$gammaPropMean <- matrix(colMeans(gammaPropSave), p, m)
  colnames(samplerPar$gammaPropMean) <- colnames_Y
  rownames(samplerPar$gammaPropMean) <- colnames_X
  colnames(samplerPar$pval) <- colnames_Y
  rownames(samplerPar$pval) <- colnames_X
  colnames(samplerPar$mlogpval) <- colnames_Y
  rownames(samplerPar$mlogpval) <- colnames_X

  ########## Posterior MCMC ##########

  digits <- 3
  postMean <- NULL
  postMean$gamma <- .display_round(matrix(colMeans(abs(BSave) > 0), p, m), digits)
  postMean$theta <- .display_round(matrix(colMeans(BSave), p, m), digits)
  postMean$G <- .display_round(apply(GSave, c(1, 2), mean), digits)
  postMean$R <- .display_round(apply(RSave, c(1, 2), mean), digits)
  postMean$D <- .display_round(colMeans(DSave), digits)
  if (hyperPar$tauRandom) {
    postMean$tau <- .display_round(colMeans(tauSave), digits)
    postMean$eta <- .display_round(mean(etaSave), digits)
  }

  colnames(postMean$gamma) <- colnames_Y; rownames(postMean$gamma) <- colnames_X
  colnames(postMean$theta) <- colnames_Y; rownames(postMean$theta) <- colnames_X
  colnames(postMean$G) <- colnames_Y; rownames(postMean$G) <- colnames_Y
  colnames(postMean$R) <- colnames_Y; rownames(postMean$R) <- colnames_Y
  names(postMean$D) <- colnames_Y

  opt <- list(seed = seed)

  endTime <- Sys.time()
  timeMR2 <- paste0(.display_round(as.numeric(difftime(endTime, startTime, units = "min")), digits = 3), "m")

  list(theta = BSave, G = GSave, R = RSave, D = DSave, tau = tauSave, eta = etaSave,
      postMean = postMean, timeMR2 = timeMR2,
      hyperPar = hyperPar, samplerPar = samplerPar, opt = opt)
}


########## BVS (Gaussian-only) ##########

BVS <- function(Y, NA_Y, X, XB, gammaCurr, BCurr, ZCurr, RCurr, invRCurr, d, hyperPar, samplerPar)
{
  mlogpval <- samplerPar$mlogpval
  aPi <- hyperPar$aPi
  bPi <- hyperPar$bPi
  tau <- hyperPar$tau
  tauRandom <- isTRUE(hyperPar$tauRandom)
  eta <- NA_real_

  n <- nrow(Y); m <- ncol(Y); p <- ncol(X)
  if (tauRandom) eta <- rgamma(1, hyperPar$tauG + m * hyperPar$tauAlpha, rate = hyperPar$tauH + sum(1 / tau))
  pStar <- p
  ind <- seq_len(m)

  gammaProp <- rep(0, length(gammaCurr))
  ZTilde <- ZCurr + XB

  for (k in seq_len(m))
  {
    NAInd <- NA_Y[[k]]
    start <- (k - 1) * p + 1
    end   <- k * p
    gammaInd <- which(gammaCurr[start:end] > 0)

    ########## Updating the prior variance of the direct effects ##########

    if (tauRandom)
      tau[k] <- 1 / rgamma(1, hyperPar$tauAlpha + 0.5 * length(gammaInd),
                           rate = eta + 0.5 * sum(BCurr[start:end]^2))
    tauK <- tau[k]

    ########## Proposing gamma ##########

    gammaPropList <- Sample_gamma(gammaCurr[start:end], pStar, mlogpval[, k], samplerPar)
    loggammaPropListRatio <- gammaPropList[[2]]
    gammaStar <- gammaPropList[[1]]

    if (length(gammaStar) > 0) {
      setOne  <- gammaStar + (k - 1) * p
      setZero <- (start:end)[-gammaStar]
    } else {
      setOne  <- NULL
      setZero <- start:end
    }

    if (!is.na(bPi)) {
      logpriorgammaRatio <- lbeta(length(gammaStar) + aPi, pStar - length(gammaStar) + bPi) -
                            lbeta(sum(gammaCurr[start:end]) + aPi, pStar - sum(gammaCurr[start:end]) + bPi)
    } else {
      logpriorgammaRatio <- dbinom(length(gammaStar), pStar, aPi, log = TRUE) -
                            dbinom(sum(gammaCurr[start:end]), pStar, aPi, log = TRUE)
    }

    if (length(gammaInd) > 0) {
      Xgamma  <- as.matrix(X[, gammaInd])
      tXgamma <- t(Xgamma)
    } else {
      Xgamma  <- matrix(0, n, 1)
      tXgamma <- t(Xgamma)
    }

    if (length(gammaStar) > 0) {
      XgammaStar  <- as.matrix(X[, gammaStar])
      tXgammaStar <- t(XgammaStar)
    } else {
      XgammaStar  <- matrix(0, n, 1)
      tXgammaStar <- t(XgammaStar)
    }

    ########## Proposing Beta (Gaussian conjugate update) ##########

    RR     <- RCurr[k, -k, drop = FALSE] %*% solve(RCurr[-k, -k])
    condSD <- sqrt(max(RCurr[k, k] - RR %*% RCurr[-k, k], .Machine$double.eps))
    condMu <- drop(RR %*% t(ZCurr[, -k, drop = FALSE]))   # length-n; precomputed once per k (efficiency: does not depend on i)

    ZTilde[, k] <- d[k] * ZCurr[, k] + XB[, k]
    const <- invRCurr[k, k] / (d[k]^2)
    T_k   <- rep(0, n)
    for (kk in ind[-k]) T_k <- T_k + ZCurr[, kk] * invRCurr[k, kk] / d[k]

    if (length(gammaStar) > 1) {
      Q    <- diag(rep(1 / tauK, length(gammaStar))) + const * crossprod(XgammaStar)
      invQ <- solve(Q)
    } else {
      Q    <- as.matrix(1 / tauK) + const * crossprod(XgammaStar)
      invQ <- solve(Q)
    }
    muProp <- invQ %*% tXgammaStar %*% (const * ZTilde[, k] + T_k)

    if (length(gammaInd) > 1) {
      Q_OLD    <- diag(rep(1 / tauK, length(gammaInd))) + const * crossprod(Xgamma)
      invQ_OLD <- solve(Q_OLD)
    } else {
      Q_OLD    <- as.matrix(1 / tauK) + const * crossprod(Xgamma)
      invQ_OLD <- solve(Q_OLD)
    }
    muProp_OLD <- invQ_OLD %*% tXgamma %*% (const * ZTilde[, k] + T_k)
    BProp <- mvrnorm(1, muProp, invQ)

    XBProp <- if (length(gammaStar) > 0) (as.matrix(X[, gammaStar]) %*% BProp)[, 1] else rep(0, n)

    ########## Updating Beta and gamma ##########

    gammaProp[setOne]  <- 1
    gammaProp[setZero] <- 0

    likeHolmes     <- -0.5 * determinant(Q_OLD, logarithm = TRUE)$modulus +
                       0.5 * t(muProp_OLD) %*% Q_OLD %*% as.matrix(muProp_OLD) -
                       0.5 * length(gammaInd) * log(tauK)
    likeHolmesProp <- -0.5 * determinant(Q, logarithm = TRUE)$modulus +
                       0.5 * t(muProp) %*% Q %*% as.matrix(muProp) -
                       0.5 * length(gammaStar) * log(tauK)
    logAccRatio <- likeHolmesProp - likeHolmes + loggammaPropListRatio + logpriorgammaRatio

    if (log(runif(1)) < logAccRatio)
    {
      if (length(gammaStar) > 0) {
        BCurr[setOne] <- BProp
        gammaCurr[setOne] <- 1
      }
      gammaCurr[setZero] <- 0
      BCurr[setZero] <- 0
      XB[, k] <- XBProp
    }

    ########## Updating Z latent variable (vectorised: conditional mean/sd don't depend on i) ##########

    obsInd  <- NAInd
    missInd <- setdiff(seq_len(n), obsInd)

    if (length(obsInd) > 0)  ZCurr[obsInd,  k] <- (Y[obsInd, k] - XB[obsInd, k]) / d[k]
    if (length(missInd) > 0) ZCurr[missInd, k] <- rnorm(length(missInd), condMu[missInd], condSD)
  }

  list(gammaProp, gammaCurr, BCurr, ZCurr, XB, tau, eta)
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
      if (runif(1) < probAdd) {
        change_j <- Qrnd(p, pgam, probMix, geomMean)
        add_j <- YXcorr_zero_rank[change_j]
        setOneProp <- union(setOne, add_j)
        pgammaProp <- length(setOneProp)
        logPropgamma <- (log(1 - probAdd) - log(probAdd)) + (-log(pgammaProp) - log(Qpdf(change_j, p, pgam, probMix, geomMean)))
      } else {
        change_j <- sample.int(pgam, 1)
        del_j <- YXcorr_rank[change_j]
        setOneProp <- setdiff(setOne, del_j)
        pgammaProp <- length(setOneProp)
        logPropgamma <- (log(probAdd) - log(1 - probAdd)) + (log(pgam) + log(Qpdf(change_j, p, pgammaProp, probMix, geomMean)))
      }
    }
    if (pgam == 0 & pgam < p) {
      change_j <- Qrnd(p, 0, probMix, geomMean)
      add_j <- YXcorr_zero_rank[change_j]
      setOneProp <- union(setOne, add_j)
      pgammaProp <- length(setOneProp)
      logPropgamma <- (log(1 - probAdd) - log(1)) + (-log(pgammaProp) - log(Qpdf(change_j, p, pgam, probMix, geomMean)))
    }
    if (pgam > 0 & pgam == p) {
      change_j <- sample.int(pgam, 1)
      del_j <- YXcorr_rank[change_j]
      setOneProp <- setdiff(setOne, del_j)
      pgammaProp <- length(setOneProp)
      logPropgamma <- (log(probAdd) - log(1)) + (log(pgam) + log(Qpdf(change_j, p, pgammaProp, probMix, geomMean)))
    }
  } else {
    change_j <- Qrnd(p, pgam, probMix, geomMean)
    add_j <- YXcorr_zero_rank[change_j]
    setOneProp <- union(setOne, add_j)
    pgammaProp <- length(setOneProp)
    change_j_OLD <- change_j
    pgammaProp_OLD <- pgammaProp
    cond <- 0
    while (cond == 0) {
      change_j <- sample.int(pgammaProp, 1)
      if (abs(change_j - change_j_OLD) > 0) cond <- 1
    }
    del_j <- YXcorr_rank[change_j]
    setOneProp <- setdiff(setOne, del_j)
    pgammaProp <- length(setOneProp)
    logPropgamma <- (-log(pgammaProp_OLD) - log(Qpdf(change_j_OLD, p, pgam, probMix, geomMean))) +
                    (log(Qpdf(change_j, p, pgammaProp, probMix, geomMean)) + log(pgammaProp_OLD))
  }
  setOneProp <- sort(setOneProp)

  list(setOneProp, logPropgamma)
}

########## Qrnd ##########

Qrnd <- function(p, pgam, u, geomMean)
{
  if (runif(1) < u) sample.int(p - pgam, 1)
  else qgeom(runif(1) * pgeom(p - pgam, 1 / geomMean), 1 / geomMean)
}

########## Qpdf ##########

Qpdf <- function(r, p, pgam, u, geomMean)
{
  u * (1 / (p - pgam)) + (1 - u) * (dgeom(r, 1 / geomMean) / pgeom(p - pgam, 1 / geomMean))
}


########## .Sample_G ##########

#' @title Sample the graphical model for the residual correlation between responses
#' @description Prefers the compiled Rcpp graph-MCMC engine (\code{generate_graph_cpp()}, see
#' \code{src/MR2_graph.cpp}) when available -- profiling shows this loop (chordality checks,
#' running-intersection cliques, junction tree, edge proposal/acceptance) accounts for the large
#' majority of MR2's runtime, since it is called \code{graphParNiter} times per MCMC iteration
#' and is dominated by many small function calls rather than heavy linear algebra: exactly the
#' regime where R's per-call interpreter overhead dominates and Rcpp gives the largest win. Falls
#' back to the pure-R \code{generate_graph_zo()} (identical Metropolis-Hastings target
#' distribution, validated against the Rcpp version on random chordal graphs; only the RNG stream
#' differs) if the compiled code is not available (for example when the R files are sourced directly)
#' @keywords internal
.Sample_G <- function(W, GCurr, samplerPar, delta = 2)
{
  burnin <- samplerPar$GPar$burnin
  niter  <- samplerPar$GPar$niter
  J   <- ncol(W)
  Phi <- diag(rep(1, J))
  WW  <- crossprod(W)

  if (exists("generate_graph_cpp", mode = "function") && is.loaded("_MR2_generate_graph_cpp"))
  {
    generate_graph_cpp(GCurr, delta, Phi, nrow(W), WW, burnin, niter)$g_next
  } else {
    generate_graph_zo(GCurr, delta, Phi, nrow(W), WW, burnin, niter)[[1]]
  }
}


########## ripcliques_to_jtree_zo ##########

ripcliques_to_jtree_zo <- function(cliques)
{
  p <- dim(cliques)[1]
  clique_sizes <- colSums(cliques)
  num_cliques <- length(which(clique_sizes != 0))
  score <- rep(0, p)
  jtree <- matrix(0, p, p)
  if (num_cliques >= 2)
  {
    for (i in 2:num_cliques) {
      for (k in 1:(i - 1)) score[k] <- sum(intersect_zo(cliques[, i], cliques[, k]))
      if (max(score) != 0) {
        j <- which(score == max(score))[1]
        jtree[i, j] <- 1; jtree[j, i] <- 1
      }
    }
  }
  jtree
}

########## sepsize_zo ##########

sepsize_zo <- function(cliques, jtree)
{
  p <- dim(cliques)[1]
  sepsize <- matrix(0, p, p)
  clique_sizes <- colSums(cliques)
  num_cliques <- length(which(clique_sizes != 0))
  if (num_cliques >= 1)
  {
    for (i in 1:num_cliques) {
      if ((i + 1) <= num_cliques) {
        for (j in (i + 1):num_cliques) {
          if (jtree[i, j] == 1) {
            sepsize[i, j] <- sum(intersect_zo(cliques[, i], cliques[, j]))
            sepsize[j, i] <- sepsize[i, j]
          }
        }
      }
    }
  }
  sepsize
}

########## reachability_graph ##########

reachability_graph <- function(g)
{
  p <- dim(g)[1]
  A <- g
  reach_graph <- matrix(0, p, p)
  for (i in 1:(p - 1)) { reach_graph <- reach_graph + A; A <- A %*% g }
  reach_graph[reach_graph > 0] <- 1
  reach_graph
}

########## next_edge_candidate ##########

next_edge_candidate <- function(g_current)
{
  p <- max(dim(g_current)[1], dim(g_current)[2])
  iota <- sample(p)
  c(iota[1], iota[2])
}

########## check_edge_delete_zo ##########

check_edge_delete_zo <- function(i, j, cliques)
{
  p <- dim(cliques)[1]
  delete_ok <- 1
  C <- rep(0, p)
  count_cliques <- 0
  clique_sizes <- colSums(cliques)
  num_cliques <- length(which(clique_sizes != 0))
  edge_col_vec <- rep(0, p); edge_col_vec[i] <- 1; edge_col_vec[j] <- 1
  if (num_cliques >= 1)
  {
    for (k in 1:num_cliques) {
      clique_k <- cliques[, k]
      if (sum(intersect_zo(edge_col_vec, clique_k)) == 2) {
        C <- clique_k
        count_cliques <- count_cliques + 1
        if (count_cliques > 1) { delete_ok <- 0; C <- rep(9, 3); break }
      }
    }
  }
  list(delete_ok, C)
}

########## find_clique_containing_zo ##########

find_clique_containing_zo <- function(a, cliques) min(which(cliques[a, ] != 0))

########## neighbours_vertex_zo ##########

neighbours_vertex_zo <- function(g, v_i) g[, v_i]

########## parents_vertex_zo ##########

parents_vertex_zo <- function(g, i)
{
  p <- dim(g)[1]
  ps <- neighbours_vertex_zo(g, i)
  ps[i:p] <- 0
  ps
}

########## check_edge_add_same_component_zo ##########

check_edge_add_same_component_zo <- function(a, b, jtree, sepsize, cliques)
{
  edge_ok <- 0; quit <- 0; fork <- 0
  CASE_same_branch <- 0; CASE_sats_on_b2_branch <- 0; locate_a2 <- 0
  p <- dim(cliques)[1]
  ab_edge_vec <- rep(0, p); ab_edge_vec[a] <- 1; ab_edge_vec[b] <- 1
  index_b <- find_clique_containing_zo(b, cliques)
  index_a <- find_clique_containing_zo(a, cliques)
  index_b2 <- max(c(index_a, index_b))
  if (index_b2 == index_b) { b2 <- b; a2 <- a } else { b2 <- a; a2 <- b }
  index <- index_b2
  next_index <- which(parents_vertex_zo(jtree, index) != 0)
  if ((index_b2 - 1) >= 1)
  {
    for (dummy in 1:(index_b2 - 1)) {
      if (sum(next_index) == 0 & locate_a2 == 0) { CASE_same_branch <- 0; break }
      clique_next_index <- cliques[, next_index]
      if (clique_next_index[a2] == 1 & sum(next_index) != 0) {
        index_a2 <- next_index; locate_a2 <- 1; CASE_same_branch <- 1; break
      }
      index <- next_index
      next_index <- which(parents_vertex_zo(jtree, index) != 0)
    }
  }
  if (CASE_same_branch == 0) index_a2 <- find_clique_containing_zo(a2, cliques)
  clique_a2_int_clique_b2 <- intersect_zo(cliques[, index_a2], cliques[, index_b2])
  s <- sum(clique_a2_int_clique_b2)
  if (s == 0) { edge_ok <- 0; C <- rep(0, p); quit <- 1 }
  if (jtree[index_a2, index_b2] == 1) { edge_ok <- 1; C <- union_zo(ab_edge_vec, clique_a2_int_clique_b2); quit <- 1 }
  if (quit == 0)
  {
    if (CASE_same_branch == 1) {
      bottom <- index_b2
      next_parent_b2 <- which(parents_vertex_zo(jtree, index_b2) != 0)
      while ((sum(next_parent_b2) > 0) & (next_parent_b2 >= index_a2)) {
        if (sepsize[next_parent_b2, bottom] == s) {
          edge_ok <- 1; C <- union_zo(ab_edge_vec, clique_a2_int_clique_b2); quit <- 1; break
        }
        bottom <- next_parent_b2
        next_parent_b2 <- which(parents_vertex_zo(jtree, bottom) != 0)
        if (length(next_parent_b2) == 0) break
      }
    } else {
      ancestors_a2 <- rep(0, index_a2); ancestors_b2 <- rep(0, index_b2)
      ancestors_a2_col_vec <- rep(0, p); ancestors_b2_col_vec <- rep(0, p)
      next_ancestor_a2 <- index_a2; next_ancestor_b2 <- index_b2
      if (index_a2 >= 1) {
        for (count in 1:index_a2) {
          ancestors_a2[count] <- next_ancestor_a2
          ancestors_a2_col_vec[next_ancestor_a2] <- 1
          next_ancestor_a2 <- which(parents_vertex_zo(jtree, next_ancestor_a2) != 0)
          if (sum(next_ancestor_a2) == 0) break
        }
      }
      if (index_b2 >= 1) {
        for (count in 1:index_b2) {
          ancestors_b2[count] <- next_ancestor_b2
          ancestors_b2_col_vec[next_ancestor_b2] <- 1
          next_ancestor_b2 <- which(parents_vertex_zo(jtree, next_ancestor_b2) != 0)
          if (sum(next_ancestor_b2) == 0) break
        }
      }
      fork_set <- intersect_zo(ancestors_a2_col_vec, ancestors_b2_col_vec)
      fork <- max(which(fork_set != 0))
      num_seps_branch_b2 <- length(which(ancestors_b2 > fork))
      if (num_seps_branch_b2 > 0) {
        for (count in 1:num_seps_branch_b2) {
          index_row <- ancestors_b2[count]; index_col <- ancestors_b2[count + 1]
          if (sepsize[index_row, index_col] == s) {
            edge_ok <- 1; C <- union_zo(ab_edge_vec, clique_a2_int_clique_b2)
            quit <- 1; CASE_sats_on_b2_branch <- 1; break
          }
        }
      }
      if (quit == 0 & CASE_sats_on_b2_branch == 0) {
        num_seps_branch_a2 <- length(which(ancestors_a2 > fork))
        if (num_seps_branch_a2 > 0) {
          for (count in 1:num_seps_branch_a2) {
            index_row <- ancestors_a2[count + 1]; index_col <- ancestors_a2[count]
            if (sepsize[index_row, index_col] == s) {
              edge_ok <- 1; C <- union_zo(ab_edge_vec, clique_a2_int_clique_b2); quit <- 1; break
            }
          }
        }
      }
    }
  }
  if (quit == 0) { edge_ok <- 0; C <- NULL }
  list(edge_ok, C)
}

########## next_graph_candidate_zo ##########

next_graph_candidate_zo <- function(g_current, jtree, sepsize, cliques, reach_graph)
{
  g_proposal <- g_current
  CASE_add <- 0; CASE_delete <- 0
  a <- b <- C <- NULL
  while (identical(g_proposal, g_current))
  {
    edge <- next_edge_candidate(g_current)
    i <- edge[1]; j <- edge[2]
    if (g_current[i, j] == 1) {
      ss <- check_edge_delete_zo(i, j, cliques)
      CASE_delete <- ss[[1]]
      if (CASE_delete == 1) { g_proposal[i, j] <- 0; g_proposal[j, i] <- 0; a <- i; b <- j; C <- ss[[2]] }
    } else {
      if (reach_graph[i, j] == 0) {
        CASE_add <- 1; a <- i; b <- j; C <- c(a, b)
        g_proposal[i, j] <- 1; g_proposal[j, i] <- 1
      } else {
        ss <- check_edge_add_same_component_zo(i, j, jtree, sepsize, cliques)
        CASE_add <- ss[[1]]
        if (CASE_add == 1) { a <- i; b <- j; C <- ss[[2]]; g_proposal[i, j] <- 1; g_proposal[j, i] <- 1 }
      }
    }
  }
  list(g_proposal, a, b, C, CASE_add, CASE_delete)
}

########## ln_h_ratio_zo ##########

ln_h_ratio_zo <- function(C, a, b, delta, Phi)
{
  C_nodes <- which(C != 0)
  D <- c(a, b)
  Sq2 <- setdiff(C_nodes, D)
  numSq2 <- length(Sq2)
  delta_star <- delta - 1
  Phi_CC <- rbind(cbind(matrix(Phi[Sq2, Sq2], length(Sq2), length(Sq2)), matrix(Phi[Sq2, D], length(Sq2), length(D))),
                  cbind(matrix(Phi[D, Sq2], length(D), length(Sq2)), matrix(Phi[D, D], length(D), length(D))))
  L <- t(chol(Phi_CC))
  L_DD <- L[c(numSq2 + 1, numSq2 + 2), c(numSq2 + 1, numSq2 + 2)]
  l_aa <- L_DD[1, 1]; l_bb <- L_DD[2, 2]
  ln_top <- (delta_star + numSq2 + 2) * (log(l_aa) + log(l_bb))
  ln_bottom <- (delta_star + numSq2 + 1) / 2 * (log(l_aa^2 * L_DD[2, 1]^2 + l_aa^2 * l_bb^2))
  ln_gam_bit <- -log(2 * sqrt(pi)) + lgamma((delta_star + numSq2 + 1) / 2) - lgamma((delta_star + numSq2 + 2) / 2)
  ln_top - ln_bottom + ln_gam_bit
}

########## generate_graph_zo ##########

generate_graph_zo <- function(g_current, delta, Phi, N, S_N, warmup, iters_to_keep)
{
  p <- max(dim(Phi)[1], dim(Phi)[2])
  iter <- iters_to_keep
  if (iter <= warmup) stop("not enough iterations specified")
  accept_count_sample <- 0
  g_next <- g_current

  for (i in seq_len(iter))
  {
    g_current <- g_next
    ss <- check_chordal(g_current)
    order <- ss[[2]]
    cliques <- chordal_to_ripcliques_zo(g_current, order)
    jtree <- ripcliques_to_jtree_zo(cliques)
    sepsize <- sepsize_zo(cliques, jtree)
    reach_graph <- reachability_graph(g_current)
    u <- runif(1)

    output <- next_graph_candidate_zo(g_current, jtree, sepsize, cliques, reach_graph)
    g_proposal <- output[[1]]; a <- output[[2]]; b <- output[[3]]; C <- output[[4]]
    CASE_add <- output[[5]]; CASE_delete <- output[[6]]

    likelihood_ratio <- 0
    if (CASE_add == 1) {
      ln_likelihood_ratio <- ln_h_ratio_zo(C, a, b, delta, Phi) - ln_h_ratio_zo(C, a, b, delta + N, Phi + S_N)
      likelihood_ratio <- exp(ln_likelihood_ratio)
    } else if (CASE_delete == 1) {
      ln_likelihood_ratio <- ln_h_ratio_zo(C, a, b, delta + N, Phi + S_N) - ln_h_ratio_zo(C, a, b, delta, Phi)
      likelihood_ratio <- exp(ln_likelihood_ratio)
    }

    accept_prob <- min(1, likelihood_ratio)
    if (u <= accept_prob) {
      g_next <- g_proposal
      if (i > warmup) accept_count_sample <- accept_count_sample + 1
    } else {
      g_next <- g_current
    }
  }

  list(g_next, accept_count_sample / iters_to_keep)
}


########## Sim_HIW ##########

Sim_HIW <- function(G, S, n)
{
  TT <- max(dim(G)[1], dim(G)[2])
  cc <- check_chordal(G)
  if (cc[[1]] == 0) stop("graph not decomp")
  order <- cc[[2]]
  cliques <- chordal_to_ripcliques_zo(G, order)
  Tost_Thi_hat <- g_constrain_zo(S, order, cliques)[[1]]
  sigma_identity <- generate_HIW_g_delta_identity_zo(G, cliques, n)[[1]]
  mysigma <- transform_g_conditional_HIW_no_mcs_zo(sigma_identity, G, cliques, n, diag(rep(1, TT)), Tost_Thi_hat)
  list(mysigma[[1]], mysigma[[2]])
}

########## setdiff_zo ##########

setdiff_zo <- function(a1, a2)
{
  b <- a1 - a2
  if (length(b) > 0) for (j in seq_along(b)) if (b[j] < 0) b[j] <- a1[j]
  b
}

########## setdiag ##########

setdiag <- function(M, v) { diag(M) <- v; M }

########## neighbours_vertex_cell ##########

neighbours_vertex_cell <- function(G, i) which(G[i, ] != 0)

########## union_zo ##########

union_zo <- function(v1, v2)
{
  b <- v1 + v2
  if (length(b) > 0) for (j in seq_along(b)) if (b[j] > 1) b[j] <- 1
  b
}

########## intersect_zo ##########

intersect_zo <- function(a1, a2) a1 * a2

########## check_chordal ##########

check_chordal <- function(g)
{
  g <- setdiag(g, 1)
  p <- dim(g)[1]
  order <- rep(0, p)
  chordal <- 1
  numbered <- 1
  order[1] <- 1
  for (i in 2:p)
  {
    capU <- setdiff(1:p, numbered)
    score <- rep(0, length(capU))
    if (length(capU) > 0)
      for (u_i in seq_along(capU)) score[u_i] <- length(intersect(neighbours_vertex_cell(g, capU[u_i]), numbered))
    u <- capU[which(score == max(score))[1]]
    numbered <- c(numbered, u)
    order[i] <- u
    pa <- intersect(neighbours_vertex_cell(g, u), order[1:(i - 1)])
    if (!identical(matrix(g[pa, pa], length(pa), length(pa)), matrix(1, length(pa), length(pa)))) { chordal <- 0; break }
  }
  list(chordal, order)
}

########## chordal_to_ripcliques_zo ##########

chordal_to_ripcliques_zo <- function(G, order)
{
  p <- dim(G)[1]
  pa <- matrix(0, p, p)
  num_pa <- rep(0, p)
  for (i in 2:p) {
    v <- order[i]
    pre_v <- rep(0, p); pre_v[order[1:(i - 1)]] <- 1
    pa_i <- intersect_zo(neighbours_vertex_zo(G, v), pre_v)
    num_pa[i] <- length(which(pa_i != 0))
    pa[, v] <- pa_i
  }
  ladder <- rep(0, p)
  for (i in 1:p) if ((i == p) | (sum(pa[, order[i]]) >= sum(pa[, order[i + 1]]))) ladder[i] <- order[i]
  cliques <- matrix(0, p, p)
  idx <- 0
  for (i in 1:p) {
    if (ladder[i] != 0) {
      v_col_i <- rep(0, p); v_col_i[ladder[i]] <- 1
      idx <- idx + 1
      cliques[, idx] <- union_zo(v_col_i, pa[, ladder[i]])
    }
  }
  cliques
}

########## seps_resids_hists_zo ##########

seps_resids_hists_zo <- function(cliques)
{
  p <- dim(cliques)[1]
  num_cliques <- max(which(colSums(cliques) != 0))
  seps <- matrix(0, p, p)
  resids <- hists <- matrix(0, p, p)
  hists[, 1] <- cliques[, 1]
  if (num_cliques >= 2) {
    for (index in 2:num_cliques) {
      hists[, index] <- union_zo(cliques[, index], hists[, index - 1])
      seps[, index] <- intersect_zo(cliques[, index], hists[, index - 1])
      resids[, index] <- setdiff_zo(cliques[, index], hists[, index - 1])
    }
  }
  list(seps, resids, hists)
}

########## generate_HIW_g_delta_identity_zo ##########

generate_HIW_g_delta_identity_zo <- function(g, cliques, delta)
{
  num_cliques <- max(which(colSums(cliques) != 0))
  ss <- seps_resids_hists_zo(cliques)
  residuals <- ss[[2]]
  num_Rj <- rep(0, num_cliques)
  p <- max(dim(g)[1], dim(g)[2])
  perfect_order <- rep(0, p)
  c_1 <- which(cliques[, 1] != 0)
  perfect_order[1:length(c_1)] <- c_1
  index_finish <- length(c_1)
  if (num_cliques >= 2) {
    for (j in 2:num_cliques) {
      Rj <- which(residuals[, j] != 0)
      num_Rj[j] <- length(Rj)
      index_start <- index_finish + 1
      index_finish <- index_start + num_Rj[j] - 1
      perfect_order[index_start:index_finish] <- Rj
    }
  }
  rev_perf <- rev(perfect_order)
  g_rev_perf <- g[rev_perf, rev_perf]
  Psi <- matrix(0, p, p)
  for (i in 1:p) {
    if (i == p) {
      Psi[p, p] <- rchisq(1, df = delta)^0.5
    } else {
      nu_i <- sum(g_rev_perf[i, (i + 1):p])
      Psi[i, i] <- rchisq(1, df = delta + nu_i)^0.5
    }
  }
  for (i in 1:(p - 1)) for (j in (i + 1):p) Psi[i, j] <- if (g_rev_perf[i, j] == 0) 0 else rnorm(1)

  K_rev_perf <- t(Psi) %*% Psi
  sigma_id_rev_perf <- solve(K_rev_perf)
  inverse_permute <- match(seq_len(p), rev_perf)
  K_id <- K_rev_perf[inverse_permute, inverse_permute]
  sigma_id <- sigma_id_rev_perf[inverse_permute, inverse_permute]
  list(sigma_id, K_id, sigma_id_rev_perf, K_rev_perf)
}

########## g_constrain_zo ##########

g_constrain_zo <- function(B, order, cliques)
{
  p <- dim(B)[1]
  seps <- seps_resids_hists_zo(cliques)[[1]]
  num_cliques <- max(which(colSums(cliques) != 0))
  sum_big_K_cliques <- sum_big_K_seps <- matrix(0, p, p)
  if (num_cliques >= 1) {
    for (j in 1:num_cliques) {
      cj <- which(cliques[, j] != 0)
      big_Kj <- matrix(0, p, p)
      big_Kj[cj, cj] <- solve(B[cj, cj])
      sum_big_K_cliques <- sum_big_K_cliques + big_Kj
    }
    for (j in 1:num_cliques) {
      sj <- which(seps[, j] != 0)
      big_Kj <- matrix(0, p, p)
      if (length(sj) > 0) big_Kj[sj, sj] <- solve(B[sj, sj])
      sum_big_K_seps <- sum_big_K_seps + big_Kj
    }
  }
  K_hat <- sum_big_K_cliques - sum_big_K_seps
  B_hat <- solve(K_hat + diag(rep(1e-8, dim(K_hat)[2])))
  list(B_hat, K_hat)
}

########## transform_g_conditional_HIW_no_mcs_zo ##########

transform_g_conditional_HIW_no_mcs_zo <- function(sigma_B, g, cliques, delta, B, D)
{
  p <- dim(g)[1]
  num_cliques <- max(which(colSums(cliques) != 0))
  ss <- seps_resids_hists_zo(cliques)
  seps <- ss[[1]]; residuals <- ss[[2]]
  num_Cj <- num_Sj <- num_Rj <- rep(0, num_cliques)
  num_Cj[1] <- sum(cliques[, 1])
  perfect_order <- rep(0, p)
  perfect_order[1:num_Cj[1]] <- which(cliques[, 1] != 0)
  index_finish <- num_Cj[1]
  if (num_cliques >= 2) {
    for (j in 2:num_cliques) {
      Rj <- which(residuals[, j] != 0)
      Sj <- which(seps[, j] != 0)
      num_Cj[j] <- length(which(cliques[, j] != 0))
      num_Rj[j] <- length(Rj); num_Sj[j] <- length(Sj)
      index_start <- index_finish + 1
      index_finish <- index_start + num_Rj[j] - 1
      perfect_order[index_start:index_finish] <- Rj
    }
  }
  rev_perf <- rev(perfect_order)
  g_rev_perf <- g[rev_perf, rev_perf]
  B_rev_perf <- B[rev_perf, rev_perf]
  D_rev_perf <- D[rev_perf, rev_perf]
  sigma_B_rev_perf <- sigma_B[rev_perf, rev_perf]
  choleskyK_B <- chol(solve(sigma_B_rev_perf))
  c1 <- num_Cj[1]
  indexC1_in_Upsilon_D <- (p - c1 + 1):p
  Upsilon_D <- matrix(0, p, p)
  B_1 <- B_rev_perf[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D]
  Q_1 <- chol(solve(B_1))
  D_1 <- D_rev_perf[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D]
  P_1 <- chol(solve(D_1))
  O_1 <- solve(Q_1) %*% P_1
  Upsilon_D[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D] <- choleskyK_B[indexC1_in_Upsilon_D, indexC1_in_Upsilon_D] %*% O_1

  if (num_cliques >= 2) {
    for (j in 2:num_cliques) {
      cj <- num_Cj[j]
      Rj <- which(residuals[, j] != 0); rj <- num_Rj[j]
      Sj <- which(seps[, j] != 0); sj <- num_Sj[j]
      indexRj_in_Upsilon_D <- if (rj > 0) sort(match(Rj, rev_perf)) else integer(0)
      indexSj_in_Upsilon_D <- if (sj > 0) sort(match(Sj, rev_perf)) else integer(0)
      indexCj_in_Upsilon_D <- c(indexRj_in_Upsilon_D, indexSj_in_Upsilon_D)

      B_j <- B_rev_perf[indexCj_in_Upsilon_D, indexCj_in_Upsilon_D]
      Q_j <- chol(solve(B_j))
      D_j <- D_rev_perf[indexCj_in_Upsilon_D, indexCj_in_Upsilon_D]
      P_j <- chol(solve(D_j))
      O_j <- solve(Q_j) %*% P_j

      indexRj_inOj <- if (rj > 0) 1:rj else NULL
      indexSj_inOj <- if ((rj + 1) <= cj) (rj + 1):cj else NULL

      Upsilon_D[indexRj_in_Upsilon_D, indexRj_in_Upsilon_D] <-
        as.matrix(choleskyK_B[indexRj_in_Upsilon_D, indexRj_in_Upsilon_D]) %*% as.matrix(O_j[indexRj_inOj, indexRj_inOj])
      Upsilon_D[indexRj_in_Upsilon_D, indexSj_in_Upsilon_D] <-
        as.matrix(choleskyK_B[indexRj_in_Upsilon_D, indexRj_in_Upsilon_D]) %*% O_j[indexRj_inOj, indexSj_inOj] +
        choleskyK_B[indexRj_in_Upsilon_D, indexSj_in_Upsilon_D] %*% as.matrix(O_j[indexSj_inOj, indexSj_inOj])
    }
  }
  K_D_rev_perf <- t(Upsilon_D) %*% Upsilon_D
  sigma_D_rev_perf <- solve(K_D_rev_perf)
  inverse_permute <- match(seq_len(p), rev_perf)
  K_D <- K_D_rev_perf[inverse_permute, inverse_permute]
  sigma_D <- sigma_D_rev_perf[inverse_permute, inverse_permute]
  list(sigma_D, K_D)
}


########## mod_BDgraphInit ##########

mod_BDgraphInit <- function(data, samplerPar,
                            n = NULL, not.cont = NULL, df.prior = 2, g.prior = 0.5,
                            threshold = 1e-08)
{
  burnin <- samplerPar$BDgraphPar$burnin
  niter <- samplerPar$BDgraphPar$niter
  cores <- samplerPar$BDgraphPar$cores

  if (niter < burnin) stop("Number of iterations must be higher than the burn-in")

  burnin <- floor(burnin)
  cores <- BDgraph::get_cores(cores = cores)
  list_S_n_p <- BDgraph::get_S_n_p(data = data, method = "ggm", n = n, not.cont = not.cont)
  S <- list_S_n_p$S
  n <- list_S_n_p$n
  p <- list_S_n_p$p
  colnames_data <- list_S_n_p$colnames_data
  b <- df.prior
  b_star <- b + n
  D <- diag(p)
  Ds <- D + S
  Ts <- chol(solve(Ds))
  g_prior <- BDgraph::get_g_prior(g.prior = g.prior, p = p)
  G <- BDgraph::get_g_start(g.start = "empty", g_prior = g_prior, p = p)
  G[g_prior == 1] <- 1
  G[g_prior == 0] <- 0
  G[lower.tri(G, diag = TRUE)] <- 0
  G <- G + t(G)
  K <- BDgraph::get_K_start(G = G, g.start = "empty", Ts = Ts, b_star = b_star, threshold = threshold)
  p_links <- K_hat <- matrix(0, p, p)

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

  list(p_links = p_links, last_K = last_K, last_graph = last_graph)
}


########## mod_BDgraph ##########

mod_BDgraph <- function(data, g.start, K.start, samplerPar,
                        n = NULL, not.cont = NULL, df.prior = 2, g.prior = 0.5,
                        threshold = 1e-08)
{
  burnin <- samplerPar$BDgraphPar$burnin
  niter <- samplerPar$BDgraphPar$niter
  cores <- samplerPar$BDgraphPar$cores

  if (niter < burnin) stop("Number of iterations must be larger than the burn-in")

  burnin <- floor(burnin)
  cores <- BDgraph::get_cores(cores = cores)
  list_S_n_p <- BDgraph::get_S_n_p(data = data, method = "ggm", n = n, not.cont = not.cont)
  S <- list_S_n_p$S
  n <- list_S_n_p$n
  p <- list_S_n_p$p
  colnames_data <- list_S_n_p$colnames_data
  b <- df.prior
  b_star <- b + n
  D <- diag(p)
  Ds <- D + S
  Ts <- chol(solve(Ds))
  g_prior <- BDgraph::get_g_prior(g.prior = g.prior, p = p)
  G <- g.start
  G[g_prior == 1] <- 1
  G[g_prior == 0] <- 0
  G[lower.tri(G, diag = TRUE)] <- 0
  G <- G + t(G)
  K <- K.start
  p_links <- K_hat <- matrix(0, p, p)

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

  list(p_links = p_links, last_K = last_K, last_graph = last_graph)
}
