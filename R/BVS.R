# 135 characters #####################################################################################################################
#' @title Bayesian Variable Selection
#' @description Internal function that performs one Markov chain Monte Carlo iteration to sample from the posterior distribution of 
#' all unknowns involved in the Bayesian Variable Selection step. For details, see \insertCite{Alexopoulos2021;textual}{MR2} and 
#' \insertCite{Zuber2023;textual}{MR2}
#'
#' @references
#' \insertAllCited{}


########## BVS ##########

BVS <- function(Y, NA_Y, X, XB, gammaCurr, BCurr, ZCurr, RCurr, invRCurr, d, 
                responseType, pFixed, 
                nCat, cutPoint, negBinomPar, 
                hyperPar, samplerPar)
{
  mlogpval <- samplerPar$mlogpval
  
  aPi <- hyperPar$aPi
  bPi <- hyperPar$aPi
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
    
    logpriorgammaRatio <- lbeta(length(gammaStarInit) + aPi, pStar - length(gammaStarInit) + bPi) - 
                          lbeta(sum(gammaCurr[start : end]) - pFixed + aPi, pStar - sum(gammaCurr[start : end]) + bPi)
    
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
