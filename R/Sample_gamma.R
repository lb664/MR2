# 135 characters #####################################################################################################################
#' @title Sampler of the binary latent variable
#' @description Internal function for the proposal distribution of the binary latent vector for each response based on 
#' \insertCite{Guan2011;textual}{MR2}. For details, see \insertCite{Alexopoulos2021;textual}{MR2}
#'
#' @references
#' \insertAllCited{}


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
