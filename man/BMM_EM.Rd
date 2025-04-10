% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BMM_EM.R
\name{BMM_EM}
\alias{BMM_EM}
\title{Beta density finite mixture model}
\usage{
BMM_EM(
  PPI,
  K = 2,
  iter_max = 1000,
  tol = 1e-06,
  plotting = 0,
  printing = FALSE,
  digits = 3,
  cex_list = list(cex = 1, cex.axis = 1.375, cex.lab = 1.375)
)
}
\arguments{
\item{PPI}{Marginal posterior probability of inclusion (mPPI) or edge marginal posterior of inclusion (ePPI) return by MR2}

\item{K}{Number of components. Default values set at \code{2}}

\item{iter_max}{Max number of iterations of the EM algorithm. Default values set at \code{1000}}

\item{tol}{Stopping criteria tolerance of the EM algorithm. Default values set at \code{1e-6}}

\item{plotting}{\itemize{
\item{\code{0} (default) No plots }
\item{\code{1} Plot after convergence }
\item{\code{2} Plot at each iteration of the EM algorithm } }}

\item{printing}{Printing parameters value at each iteration of the EM algorithm. Default values set at \code{FALSE}}

\item{digits}{Integer indicating the number of decimal places}

\item{cex_list}{List of printing options used \code{list(cex, cex.axis, cex.lab)}}
}
\value{
The value returned is a list object \code{list(PPI, w, a, b, p, Llik, aStore, bStore, pStore, LlikStore, stopcode, opt)}
\itemize{
\item{\code{PPI}}{ Ordered vector of the posterior probabilities of inclusion }
\item{\code{w}}{ Matrix of allocation probabilities }
\item{\code{a}}{ Vector of estimated first parameter of the beta components }
\item{\code{b}}{ Vector of estimated second parameter of the beta components }
\item{\code{Llik}}{ Final log-likelihood estimated by the EM algorithm }
\item{\code{aStore}}{ Stored vectors of the estimated first parameter of the beta components for each iteration of the EM
algorithm }
\item{\code{bStore}}{ Stored vectors of the estimated second parameter of the beta components for each iteration of the EM
algorithm }
\item{\code{pStore}}{ Stored vectors of estimated probability of each beta component for each iteration of the EM algorithm }
\item{\code{LlikStoreStore}}{ Stored value of the log-likelihood for each iteration of the EM algorithm }
\item{\code{stopcode}}{ Message that describes how the EM algorithm has terminated }
\item{\code{opt}}{ List of options used \code{list(iter_max, tol, plotting, printing)} } }
}
\description{
Post-processing of MR2 Markov chain Monte Carlo output. Beta density finite mixture model of the marginal posterior
probability of inclusion (mPPI) or edge marginal posterior of inclusion (ePPI)
}
\details{
For details regarding the EM algorithm to classify posterior probabilities of inclusion, see in \insertCite{Zuber2023;textual}{MR2}
}
\examples{
# Example 1: Analysis of one replication of simulated Scenario II-Confounding with q = 5 
# responses, p = 15 exposures and nIV = 100 genetic variants used as IVs. The number of 
# expected exposures directly associated with each response is set at 2 and its variance at 2, 
# with a priori range of direct causal association ranging between 0 and 8

Sim_Confounding <- Scenario_Confounding(nRep = 1, q = 5, p = 15, nIV = 100, seed = 280610971)
betaHat_Y <- Sim_Confounding$betaHat_Y
betaHat_X <- Sim_Confounding$betaHat_X

MR2_output <- MR2(betaHat_Y, betaHat_X, EVgamma = c(2, 2), 
                  niter = 7500, burnin = 2500, thin = 5, monitor = 1000, seed = 280610971)

PostProc_output <- PostProc(MR2_output, betaHat_Y, betaHat_X)

BMM_EM_output <- BMM_EM(PostProc_output$gammaPost)


# Example 2: Analysis of one replication of simulated Scenario IV-Directed pleiotropy with q = 
# 5 responses, p = 15 exposures and nIV = 100 genetic variants used as IVs and the effect of 
# the shared pleiotropic pathway on the responses set as \code{2}. The number of expected 
# exposures directly associated with each response is set at 1 and its variance at 2, 
# with a priori range of direct causal associations ranging between 0 and 7. A non-
# decomposable graph for the inverse of the residual correlation between responses is selected

Sim_Pleiotropy <- Scenario_Pleiotropy(nRep = 1, q = 5, p = 15, nIV = 100, undirectedA = FALSE, 
                                      seed = 280610971)
betaHat_Y <- Sim_Pleiotropy$betaHat_Y
betaHat_X <- Sim_Pleiotropy$betaHat_X

MR2_output <- MR2(betaHat_Y, betaHat_X, EVgamma = c(1, 3),
                  niter = 15000, burnin = 5000, thin = 10, monitor = 1000, 
                  nonDecomp = TRUE, seed = 280610971)

PostProc_output <- PostProc(MR2_output, betaHat_Y, betaHat_X)

BMM_EM_output <- BMM_EM(PostProc_output$gammaPost)
}
\references{
\insertAllCited{}
}
