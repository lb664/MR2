# 135 characters #####################################################################################################################
#' @title Beta density finite mixture model
#' @description Post-processing of MR2 Markov chain Monte Carlo output. Beta density finite mixture model of the marginal posterior 
#' probability of inclusion (mPPI) or edge marginal posterior of inclusion (ePPI)
#'
#' @param PPI Marginal posterior probability of inclusion (mPPI) or edge marginal posterior of inclusion (ePPI) return by MR2
#' @param K Number of components. Default values set at \code{2}
#' @param iter_max Max number of iterations of the EM algorithm. Default values set at \code{1000}
#' @param tol Stopping criteria tolerance of the EM algorithm. Default values set at \code{1e-6}
#' @param plotting \itemize{
#'   \item{\code{0} (default) No plots }
#'   \item{\code{1} Plot after convergence }
#'   \item{\code{2} Plot at each iteration of the EM algorithm } }
#' @param printing Printing parameters value at each iteration of the EM algorithm. Default values set at \code{FALSE}
#' @param digits Integer indicating the number of decimal places
#' @param cex_list List of printing options used \code{list(cex, cex.axis, cex.lab)}
#'
#' @details
#' For details regarding the EM algorithm to classify posterior probabilities of inclusion, see in \insertCite{Zuber2023;textual}{MR2}
#'
#' @export
#'
#' @return The value returned is a list object \code{list(PPI, w, a, b, p, Llik, aStore, bStore, pStore, LlikStore, stopcode, opt)}
#' \itemize{
#'   \item{\code{PPI}}{ Ordered vector of the posterior probabilities of inclusion }
#'   \item{\code{w}}{ Matrix of allocation probabilities }
#'   \item{\code{a}}{ Vector of estimated first parameter of the beta components }
#'   \item{\code{b}}{ Vector of estimated second parameter of the beta components }
#'   \item{\code{Llik}}{ Final log-likelihood estimated by the EM algorithm }
#'   \item{\code{aStore}}{ Stored vectors of the estimated first parameter of the beta components for each iteration of the EM 
#'         algorithm }
#'   \item{\code{bStore}}{ Stored vectors of the estimated second parameter of the beta components for each iteration of the EM 
#'         algorithm }
#'   \item{\code{pStore}}{ Stored vectors of estimated probability of each beta component for each iteration of the EM algorithm }
#'   \item{\code{LlikStoreStore}}{ Stored value of the log-likelihood for each iteration of the EM algorithm }
#'   \item{\code{stopcode}}{ Message that describes how the EM algorithm has terminated }
#'   \item{\code{opt}}{ List of options used \code{list(iter_max, tol, plotting, printing)} } }
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
# 100 characters ##################################################################################
#' # Example 1: Analysis of one replication of simulated Scenario II-Confounding with q = 5 
#' # responses, p = 15 exposures and nIV = 100 genetic variants used as IVs. The number of 
#' # expected exposures directly associated with each response is set at 2 and its variance at 2, 
#' # with a priori range of direct causal association ranging between 0 and 8
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
#' # with a priori range of direct causal associations ranging between 0 and 7. A non-
#' # decomposable graph for the inverse of the residual correlation between responses is selected
#'
#' Sim_Pleiotropy <- Scenario_Pleiotropy(nRep = 1, q = 5, p = 15, nIV = 100, undirectedA = FALSE, 
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


BMM_EM <- function(PPI, K = 2, iter_max = 1000, tol = 1e-6, plotting = 0, printing = FALSE, digits = 3, 
                   cex_list = list(cex = 1, cex.axis = 1.375, cex.lab = 1.375))
{
    
##########  Functions ########## 
    
    f <- function(par, PPI, p, w)
    {
        N <- length(PPI)
        K <- length(par) /2
        d <- matrix(0, N, K)
        
        for (k in 1 : K)
        {
            a_tmp <- par[1 + 2 * (k - 1)]
            b_tmp <- par[2 + 2 * (k - 1)]
            d[, k] <- dbeta(PPI, a_tmp, b_tmp)
        }
        
        d[d == 0] <- 1e-6   # ***
        output <- -sum(w * ((rep(1, N) %*% log(p)) + log(d)))
        
        return(output)
    }
    
    
#######################
# Data transformation #
#######################
    
    # LL <- 1e-6   # ***
    # UL <- 1e-4   # ***
    # PPI[PPI == 0] <- seq(LL, UL, length.out = sum(PPI == 0))
    # PPI[PPI == 1] <- seq(1 - UL, 1 - LL, length.out = sum(PPI == 1))
    PPI_OLD <- sort(PPI, decreasing = FALSE)
    PPI <- (PPI * (length(PPI) - 1) + 0.5) / length(PPI)
    PPI <- sort(PPI, decreasing = FALSE)
    
    
##################
# Initial values #
##################
    
    N <- length(PPI)
    a <- rep(1, K)
    b <- rep(1, K)
    p <- rep(1 /K, K)
    
    output <- kmeans(PPI, K)
    
    for (k in 1 : K)
    {
        cl_idx <- output$cluster == k
        mu <- mean(PPI[cl_idx])
        sigma2 <- var(PPI[cl_idx])
        
        if (is.na(sigma2) | sigma2 == 0)
        {
            a[k] <- 1 /2
            b[k] <- 1 /2
        } else {
            a[k] <- mu * (mu * (1 - mu) / sigma2 - 1)
            b[k] <- (1 - mu) * (mu * (1 - mu) / sigma2 - 1)
        }
        
        p[k] <- sum(cl_idx) / N
    }
    
    
##############
# Reordering #
##############
    
    mu <- a / (a + b)
    mu_idx <- sort(mu, index.return = TRUE)$ix
    a <- a[mu_idx]
    b <- b[mu_idx]
    p <- p[mu_idx]
    
    
###########
# Storing #
###########
    
    aStore <- a
    bStore <- b
    pStore <- p
    Llik <- -1e10   # ***
    LlikStore <- Llik
    
    w <- matrix(0, N, K)
    d <- matrix(0, N, K)
    ww <- matrix(0, N, K)
    
    iter <- 1
    stopcond <- 0
    stopcode <- NA
    
    while (stopcond != 1)
    {
        w0 <- w
        a0 <- a
        b0 <- b
        p0 <- p
        Llik0 <- Llik
        
        
##########
# E-step #
##########
        
        for (k in 1 : K)
        {
            r <- dbeta(PPI, a[k], b[k])
            w[, k] <- p[k] * r
            d[, k] <- r
        }
	
	w <- w / (rowSums(w) %*% t(rep(1, k)))
        p <- t(colMeans(w))
        
        for (k in 1 : K)
        {
            ww[, k] <- p[k] * d[, k]
        }
        
        Llik <- sum(log(rowSums(ww)))
        
        
##########
# M-step #
##########
        
        par <- as.numeric(c(a, b))
        options(warn = - 1)
        nlm_output <- nlm(f, par, PPI, p, w, print.level = 0)
        options(warn = 0)
        
        for (k in 1 : K)
        {
            a[k] <- nlm_output$estimate[1 + 2 * (k - 1)]
            b[k] <- nlm_output$estimate[2 + 2 * (k - 1)]
        }
        
        
#####################
# Breaking criteria #
#####################
        
        if (any(is.na(c(a, b, p))))
        {
            w <- w0
            a <- a0
            b <- b0
            p <- p0
            Llik <- Llik0
            
            stopcode <- "Warning: NA beta parameters"
            stopcond <- 1
            cat(stopcode, "\n\n")
            
            break
        }
        
        if (max(abs(a / (a + b) - a0 / (a0 + b0))) > 1 /K)
        {
            w <- w0
            a <- a0
            b <- b0
            p <- p0
            Llik <- Llik0
            
            stopcode <- "Warning: M-step convergence"
            stopcond <- 1
            cat(stopcode, "\n\n")
            
            break
        }
        
        if ((Llik - Llik0) < 0)
        {
            w <- w0
            a <- a0
            b <- b0
            p <- p0
            Llik <- Llik0
            
            stopcode <- "Warning: Decreasing log-likelihood"
            stopcond <- 1
            cat(stopcode, "\n\n")
            
            break
        }
        
        
###########
# Storing #
###########
        
	aStore <- rbind(aStore, a)
	bStore <- rbind(bStore, b)
	pStore <- rbind(pStore, p)
        LlikStore <- c(LlikStore, Llik)
        
        
#####################
# Stopping criteria #
#####################
        
        dif_par <- max(abs(c(a, b, p) - c(a0, b0, p0)))
        dif_Llik <- abs(Llik - Llik0)
	
        if (iter > iter_max)
        {
            stopcode <- "Warning: Max number iterations reached"
            stopcond <- 1
            cat(stopcode, "\n\n")
        }
	
        if (iter <= iter_max)
        {
            
            if (dif_par < tol | dif_Llik < tol)
            {
                stopcode <- "Convergence reached"
                stopcond <- 1
                cat(stopcode, "\n\n")
            }
        
        }
        
        
############
# Plotting #
############

        if (plotting != 0)
        {
            cex <- cex_list$cex
            cex.axis <- cex_list$cex.axis
            cex.lab <- cex_list$cex.lab
        }
        
        if (plotting == 1)
        {
            cl <- apply(w, 1, which.max)
            PPI_max <- max(d)
            mu <- a / (a + b)
            col <- gray(seq(0.2, 0.8, length = 4))
            
            for (k in 1 : K)
            {
                
                if (k == 1)
                {
                    hist(PPI[cl == k], breaks = seq(0, 1, length.out = 10), freq = FALSE, 
                         col = col[k], main = "", xlim = c(0, 1), xlab = "PPI", ylim = c(0, PPI_max), axes = FALSE, 
                         cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.lab)
                    axis(1, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5), 
                         cex = cex, cex.axis = cex.axis, cex.lab = cex.lab)
                    axis(2, cex = cex, cex.axis = cex.axis, cex.lab = cex.lab)
                } else {
                    hist(PPI[cl == k], breaks = seq(0, 1, length.out = 10), freq = FALSE, 
                         col = col[k], add = TRUE)
                }
            
            }
            
            for (k in 1 : K)
            {
                PPIs <- seq(0, 1, by = 0.001)
                lines(PPIs, dbeta(PPIs, shape1 = a[k], shape2 = b[k]), 
                      col = col[k], lwd = 3)
                abline(v = mu[k], col = col[k], lty = 2, lwd = 3)
            }
        
        }
        
        
############
# Printing #
############
        
        if (printing == TRUE)
        {
            mu <- a / (a + b)
            sigma2 <- (a * b) / ((a + b) ^2 + (a + b + 1))
            mu_idx <- sort(mu, index.return = TRUE)$ix
            tmp <- rbind(p, a, b, mu, sigma2)
            tmp <- tmp[, mu_idx]
            rownames(tmp) <- c("Prob.", "a", "b", "mu", "sigma2")
            colnames(tmp) <- c(paste("Comp.", seq(1, K)))
            cat("BMM EM iteration:", iter, "\n")
            cat(round(tmp, digits), "\n")
        }
    
    iter <- iter + 1
    }
    
    
############
# Plotting #
############
    
    if (plotting == 2)
    {
        cl <- apply(w, 1, which.max)
        PPI_max <- max(d)
        mu <- a / (a + b)
        col <- gray(seq(0.2, 0.8, length = 4))
        
        for (k in 1 : K)
        {
            
            if (k == 1)
            {
                hist(PPI[cl == k], breaks = seq(0, 1, length.out = 10), freq = FALSE, 
                     col = col[k], main = "", xlim = c(0, 1), xlab = "PPI", ylim = c(0, PPI_max), 
                     axes = FALSE, cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.lab)
                axis(1, at = seq(0, 1, length.out = 5), labels = seq(0, 1, length.out = 5), 
                     cex = cex, cex.axis = cex.axis, cex.lab = cex.lab)
                     axis(2, cex = cex, cex.axis = cex.axis, cex.lab = cex.lab)
            } else {
                hist(PPI[cl == k], breaks = seq(0, 1, length.out = 10), freq = FALSE, 
                     col = col[k], add = TRUE)
            }
        
        }
        
        for (k in 1 : K)
        {
            PPIs <- seq(0, 1, by = 0.001)
            lines(PPIs, dbeta(PPIs, shape1 = a[k], shape2 = b[k]), 
                  col = col[k], lwd = 3)
            abline(v = mu[k], col = col[k], lty = 2, lwd = 3)
        }
    
    }
    
    
##############
# Reordering #
##############
    
    mu <- a / (a + b)
    mu_idx <- sort(mu, index.return = TRUE)$ix
    w <- w[, mu_idx]
    a <- a[mu_idx]
    b <- b[mu_idx]
    p <- p[mu_idx]
    
    if (is.matrix(aStore))
    {
        aStore <- aStore[, mu_idx]
        bStore <- bStore[, mu_idx]
        pStore <- pStore[, mu_idx]
    } else {
        aStore <- rbind(c(NA, NA), aStore[mu_idx])
        bStore <- rbind(c(NA, NA), bStore[mu_idx])
        pStore <- rbind(c(NA, NA), pStore[mu_idx])
    }
    
    names(aStore) <- NULL
    names(bStore) <- NULL
    names(pStore) <- NULL
    names(LlikStore) <- NULL
    
    opt <- list(iter_max = iter_max, tol = tol, plotting = plotting, printing = printing)
    
    output <- list(PPI = PPI_OLD, 
                   w = w, 
                   a = a, 
                   b = b, 
                   p = p, 
                   Llik = Llik, 
                   aStore = aStore[-1, ], 
                   bStore = bStore[-1, ], 
                   pStore = pStore[-1, ], 
                   LlikStore = LlikStore[-1], 
                   stopcode = stopcode, 
                   opt = opt)
    
    return(output)
}
