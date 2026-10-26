# 135 characters #####################################################################################################################
#' @title Beta density finite mixture model
#' @description Post-processing of MR2 Markov chain Monte Carlo output. Beta density finite mixture model of the marginal posterior
#' inclusion probability (mPIP) or posterior edge inclusion probability (PEIP)
#'
#' @param PIP The mPIP or PEIP returned by
#' \code{PostProc()}
#' @param K Number of components. Default value set at \code{2}
#' @param iter_max Maximum number of iterations of the EM algorithm. Default value set at \code{1000}
#' @param tol Stopping criterion tolerance of the EM algorithm. Default value set at \code{1e-6}
#' @param plotting \itemize{
#'   \item{\code{0} (default) No plots }
#'   \item{\code{1} Plot at each iteration of the EM algorithm }
#'   \item{\code{2} Plot after convergence } }
#' @param printing Print parameter values at each iteration of the EM algorithm. Default value set at \code{FALSE}
#' @param digits Integer indicating the number of decimal places
#' @param cex_list List of plotting options used \code{list(cex, cex.axis, cex.lab)}
#'
#' @details
#' For details regarding the Beta density finite mixture model, see \insertCite{Broet2004;textual}{MR2}. For details
#' regarding its use to classify the mPIP and PEIP in \strong{MR2}, see \insertCite{Zuber2023;textual}{MR2}
#'
#' @export
#'
#' @return The value returned is a list object \code{list(PIP, w, a, b, p, Llik, aStore, bStore, pStore, LlikStore, stopcode, opt)}
#' \describe{
#'   \item{\code{PIP}}{ Ordered vector of the mPIP or PEIP }
#'   \item{\code{w}}{ Matrix of allocation probabilities }
#'   \item{\code{a}}{ Vector of estimated first parameters of the beta components }
#'   \item{\code{b}}{ Vector of estimated second parameters of the beta components }
#'   \item{\code{p}}{ Vector of estimated mixing proportions of the beta components }
#'   \item{\code{Llik}}{ Final log-likelihood estimated by the EM algorithm }
#'   \item{\code{aStore}}{ Stored vectors of the estimated first parameters of the beta components for each iteration }
#'   \item{\code{bStore}}{ Stored vectors of the estimated second parameters of the beta components for each iteration }
#'   \item{\code{pStore}}{ Stored vectors of estimated mixing probabilities for each iteration }
#'   \item{\code{LlikStore}}{ Stored values of the log-likelihood for each iteration }
#'   \item{\code{stopcode}}{ Message describing how the EM algorithm terminated }
#'   \item{\code{opt}}{ List of options used \code{list(iter_max, tol, plotting, printing)} } }
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
#'
#'
#' # Example 2: Analysis of one replication of simulated Scenario IV-Directed pleiotropy with q =
#' # 5 responses, p = 15 exposures and nIV = 100 genetic variants used as IVs and the effect of
#' # the shared pleiotropic pathway on the responses set as \code{2}. The number of expected
#' # exposures directly associated with each response is set at 1 and its variance at 2,
#' # with a priori range of direct associations ranging between 0 and 7. A non-
#' # decomposable graph for the inverse of the residual correlation between responses is selected
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
#' }


BMM_EM <- function(PIP, K = 2, iter_max = 1000, tol = 1e-6, plotting = 0, printing = FALSE,
                   digits = 3, cex_list = list(cex = 1, cex.axis = 1.375, cex.lab = 1.375))
{

########## Internal objective for M-step ##########

  f_obj <- function(par, PIP, p, w)
  {
    N <- length(PIP)
    K <- length(par) / 2
    d <- matrix(0, N, K)

    for (k in seq_len(K))
    {
      a_tmp <- par[1 + 2 * (k - 1)]
      b_tmp <- par[2 + 2 * (k - 1)]
      d[, k] <- dbeta(PIP, a_tmp, b_tmp)
    }

    d[d == 0] <- 1e-6
    -sum(w * (matrix(log(p), N, K, byrow = TRUE) + log(d)))
  }


#######################
# Data transformation #
#######################

  PIP_OLD <- sort(PIP, decreasing = FALSE)
  PIP <- (PIP * (length(PIP) - 1) + 0.5) / length(PIP)
  PIP <- sort(PIP, decreasing = FALSE)


##################
# Initial values #
##################

  N <- length(PIP)
  a <- rep(1, K)
  b <- rep(1, K)
  p <- rep(1 / K, K)

  km_out <- kmeans(PIP, K)

  for (k in seq_len(K))
  {
    cl_idx <- km_out$cluster == k
    mu_k   <- mean(PIP[cl_idx])
    s2_k   <- var(PIP[cl_idx])

    if (is.na(s2_k) || s2_k == 0) {
      a[k] <- 0.5
      b[k] <- 0.5
    } else {
      common <- mu_k * (1 - mu_k) / s2_k - 1
      a[k]   <- mu_k * common
      b[k]   <- (1 - mu_k) * common
    }

    p[k] <- sum(cl_idx) / N
  }


##############
# Reordering #
##############

  mu      <- a / (a + b)
  mu_idx  <- order(mu)
  a <- a[mu_idx]; b <- b[mu_idx]; p <- p[mu_idx]


###########
# Storing #
###########

  aStore    <- a
  bStore    <- b
  pStore    <- p
  Llik      <- -1e10
  LlikStore <- Llik

  w   <- matrix(0, N, K)
  d   <- matrix(0, N, K)
  ww  <- matrix(0, N, K)

  iter     <- 1
  stopcond <- 0
  stopcode <- NA

  while (stopcond != 1)
  {
    w0 <- w; a0 <- a; b0 <- b; p0 <- p; Llik0 <- Llik


##########
# E-step #
##########

    for (k in seq_len(K))
    {
      r       <- dbeta(PIP, a[k], b[k])
      w[, k]  <- p[k] * r
      d[, k]  <- r
    }

    ## Row-normalise responsibility matrix
    w <- sweep(w, 1, rowSums(w), "/")
    p <- colMeans(w)

    for (k in seq_len(K))
      ww[, k] <- p[k] * d[, k]

    Llik <- sum(log(rowSums(ww)))


##########
# M-step #
##########

    par_vec <- as.numeric(c(a, b))
    options(warn = -1)
    nlm_out <- nlm(f_obj, par_vec, PIP, p, w, print.level = 0)
    options(warn = 0)

    for (k in seq_len(K))
    {
      a[k] <- nlm_out$estimate[1 + 2 * (k - 1)]
      b[k] <- nlm_out$estimate[2 + 2 * (k - 1)]
    }


#####################
# Breaking criteria #
#####################

    if (any(is.na(c(a, b, p))))
    {
      if (iter > 1) { w <- w0; Llik <- Llik0 }   # at the first iteration w0 is empty: keep the E-step output
      a <- a0; b <- b0; p <- p0
      stopcode <- "Warning: NA beta parameters"
      stopcond <- 1
      cat(stopcode, "\n\n")
      break
    }

    if (max(abs(a / (a + b) - a0 / (a0 + b0))) > 1 / K)
    {
      if (iter > 1) { w <- w0; Llik <- Llik0 }   # at the first iteration w0 is empty: keep the E-step output
      a <- a0; b <- b0; p <- p0
      stopcode <- "Warning: M-step convergence"
      stopcond <- 1
      cat(stopcode, "\n\n")
      break
    }

    if ((Llik - Llik0) < 0)
    {
      if (iter > 1) { w <- w0; Llik <- Llik0 }   # at the first iteration w0 is empty: keep the E-step output
      a <- a0; b <- b0; p <- p0
      stopcode <- "Warning: Decreasing log-likelihood"
      stopcond <- 1
      cat(stopcode, "\n\n")
      break
    }


###########
# Storing #
###########

    aStore    <- rbind(aStore, a)
    bStore    <- rbind(bStore, b)
    pStore    <- rbind(pStore, p)
    LlikStore <- c(LlikStore, Llik)


#####################
# Stopping criteria #
#####################

    dif_par  <- max(abs(c(a, b, p) - c(a0, b0, p0)))
    dif_Llik <- abs(Llik - Llik0)

    if (iter > iter_max)
    {
      stopcode <- "Warning: Max number of iterations reached"
      stopcond <- 1
      cat(stopcode, "\n\n")
    }

    if (iter <= iter_max && (dif_par < tol || dif_Llik < tol))
    {
      stopcode <- "Convergence reached"
      stopcond <- 1
      cat(stopcode, "\n\n")
    }


############
# Plotting #
############
    ## plotting == 1 : plot at every iteration

    if (plotting == 1)
    {
      cex     <- cex_list$cex
      cex.ax  <- cex_list$cex.axis
      cex.lab <- cex_list$cex.lab

      cl      <- apply(w, 1, which.max)
      PIP_max <- max(d)
      mu      <- a / (a + b)
      col     <- gray(seq(0.2, 0.8, length.out = 4))

      for (k in seq_len(K))
      {
        if (k == 1) {
          hist(PIP[cl == k], breaks = seq(0, 1, length.out = 10), freq = FALSE,
               col = col[k], main = "", xlim = c(0, 1), xlab = "PIP",
               ylim = c(0, PIP_max), axes = FALSE,
               cex = cex, cex.axis = cex.ax, cex.lab = cex.lab, cex.main = cex.lab)
          axis(1, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5),
               cex = cex, cex.axis = cex.ax, cex.lab = cex.lab)
          axis(2, cex = cex, cex.axis = cex.ax, cex.lab = cex.lab)
        } else {
          hist(PIP[cl == k], breaks = seq(0, 1, length.out = 10), freq = FALSE,
               col = col[k], add = TRUE)
        }
      }

      for (k in seq_len(K))
      {
        PIPs <- seq(0, 1, by = 0.001)
        lines(PIPs, dbeta(PIPs, shape1 = a[k], shape2 = b[k]), col = col[k], lwd = 3)
        abline(v = mu[k], col = col[k], lty = 2, lwd = 3)
      }
    }


############
# Printing #
############

    if (printing)
    {
      mu     <- a / (a + b)
      ## Correct Beta variance: ab / [(a+b)^2 * (a+b+1)]
      sigma2 <- (a * b) / ((a + b)^2 * (a + b + 1))
      mu_idx <- order(mu)
      tmp    <- rbind(p, a, b, mu, sigma2)[, mu_idx]
      rownames(tmp) <- c("Prob.", "a", "b", "mu", "sigma2")
      colnames(tmp) <- paste("Comp.", seq_len(K))
      cat("BMM EM iteration:", iter, "\n")
      print(.display_round(tmp, digits))
      cat("\n")
    }

    iter <- iter + 1
  }


############
# Plotting #
############
  ## plotting == 2 : plot after convergence

  if (plotting == 2)
  {
    cex     <- cex_list$cex
    cex.ax  <- cex_list$cex.axis
    cex.lab <- cex_list$cex.lab

    cl      <- apply(w, 1, which.max)
    PIP_max <- max(d)
    mu      <- a / (a + b)
    col     <- gray(seq(0.2, 0.8, length.out = 4))

    for (k in seq_len(K))
    {
      if (k == 1) {
        hist(PIP[cl == k], breaks = seq(0, 1, length.out = 10), freq = FALSE,
             col = col[k], main = "", xlim = c(0, 1), xlab = "PIP",
             ylim = c(0, PIP_max), axes = FALSE,
             cex = cex, cex.axis = cex.ax, cex.lab = cex.lab, cex.main = cex.lab)
        axis(1, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5),
             cex = cex, cex.axis = cex.ax, cex.lab = cex.lab)
        axis(2, cex = cex, cex.axis = cex.ax, cex.lab = cex.lab)
      } else {
        hist(PIP[cl == k], breaks = seq(0, 1, length.out = 10), freq = FALSE,
             col = col[k], add = TRUE)
      }
    }

    for (k in seq_len(K))
    {
      PIPs <- seq(0, 1, by = 0.001)
      lines(PIPs, dbeta(PIPs, shape1 = a[k], shape2 = b[k]), col = col[k], lwd = 3)
      abline(v = a[k] / (a[k] + b[k]), col = col[k], lty = 2, lwd = 3)
    }
  }


##############
# Reordering #
##############

  mu      <- a / (a + b)
  mu_idx  <- order(mu)
  w <- w[, mu_idx, drop = FALSE]
  a <- a[mu_idx]; b <- b[mu_idx]; p <- p[mu_idx]

  if (is.matrix(aStore)) {
    aStore <- aStore[, mu_idx, drop = FALSE]
    bStore <- bStore[, mu_idx, drop = FALSE]
    pStore <- pStore[, mu_idx, drop = FALSE]
  } else {
    aStore <- rbind(c(NA, NA), aStore[mu_idx])
    bStore <- rbind(c(NA, NA), bStore[mu_idx])
    pStore <- rbind(c(NA, NA), pStore[mu_idx])
  }

  names(aStore)    <- NULL
  names(bStore)    <- NULL
  names(pStore)    <- NULL
  names(LlikStore) <- NULL

  opt <- list(iter_max = iter_max, tol = tol, plotting = plotting, printing = printing)

  output <- list(PIP       = PIP_OLD,
                 w         = w,
                 a         = a,
                 b         = b,
                 p         = p,
                 Llik      = Llik,
                 aStore    = aStore[-1, , drop = FALSE],
                 bStore    = bStore[-1, , drop = FALSE],
                 pStore    = pStore[-1, , drop = FALSE],
                 LlikStore = LlikStore[-1],
                 stopcode  = stopcode,
                 opt       = opt)

  return(output)
}
