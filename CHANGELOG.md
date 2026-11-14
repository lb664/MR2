
# News/Changelog

Changes in **MR<sup>2</sup>** package version **0.2.0** (25/10/2026):
<ul>
    <li> Fixed a bug in BMM_EM(): when the EM algorithm stopped at its first iteration (stop code 'M-step convergence'), the matrix of allocation probabilities 'w' was returned empty (all zeros) and the log-likelihood as -1e10, so that FDR_BMM() estimated a zero false discovery rate and declared (almost) all exposure-response pairs significant. The allocation probabilities computed at the first iteration are now kept
    <li> In BMM_EM() and FDR_BMM(), the argument 'PPI' (and the output element '$PPI' of BMM_EM()) is renamed 'PIP', since it holds either the marginal posterior inclusion probability (mPIP) or the posterior edge inclusion probability (PEIP). Code that uses 'BMM_EM_output$PPI' needs 'BMM_EM_output$PIP'
    <li> The arguments of the four Scenario_*() simulators have been renamed (heritY -> herit_Y, heritX -> herit_X, rhoX -> rho_X, rhoY -> rho_Y, thetaUY -> psi_Y, thetaUX -> psi_X, thetaM -> gamma_Y, undirectedA -> undirected_A, thetaA -> delta_A). Code written for 0.1.2 needs the new names
    <li> MR2() again places a hyperprior on the prior variance of the direct effects (in version 0.1.2 and in the previous versions of the package the prior variance was fixed): one 'tau_k' per response with tau_k ~ Inverse-Gamma(alpha, eta) and eta ~ Gamma(g, h), updated by Gibbs sampling at every MCMC iteration (hyperparameters 'tauPar = c(alpha = 2, g = 0.2, h = 10)'). This is the new default ('tau = NULL'); a single positive number, e.g. 'tau = 1' (the previous default), still fixes the prior variance for every response. The posterior draws of 'tau' (iterations x responses) and 'eta' are returned as '$tau' and '$eta', and their posterior means are added to '$postMean'
    <li> The graphical model Markov chain Monte Carlo sampler (Sample_G(), the dominant cost of MR2()) can now be accelerated with a compiled Rcpp/RcppArmadillo implementation ('src/MR2_graph.cpp')
    <li> Added 'Rcpp' to Imports and 'Rcpp'/'RcppArmadillo' to LinkingTo
    <li> Added internal 'display_round()' (see 'R/utils.R'): rounds to 'digits' decimals as 'round()' does, but falls back to 'signif()' for elements that would otherwise display as an exact zero. Used throughout MR2()'s and PostProc()'s human-facing summary output ('postMean', 'gammaPost', 'thetaPost', 'DPost', 'GPost', 'SigmaPost', 'RPost', 'Rm1Post', 'PLML', 'DIC', timing strings, etc.) and in BMM_EM()'s iteration printout, so small non-zero estimates are no longer misread as exactly zero
    <li> MR2()'s MCMC progress is now shown as a live console progress bar (refreshed every 'monitor'-th iteration) instead of printed "iteration N" text
    <li> Fixed a bug where 'BSave'/'GSave'/'RSave'/'DSave' were sized as '(niter - burnin) / thin', silently truncated to an integer by matrix()/array() whenever that division was not exact (e.g. niter = 20, burnin = 5, thin = 2 gives 7.5 -> 7 rows allocated for 8 actual saves), causing "subscript out of bounds" partway through the MCMC run. These are now sized from the exact set of iterations satisfying 'iter > burnin & iter %% thin == 0'
    <li> Removed the 'fullCov' option from MR2(): the dense (full-covariance, non-graphical) inverse-correlation sampler has been dropped entirely, along with its 'MCMCpack::rwish()' dependency (removed from Imports). MR2() now always uses the sparse graphical model, decomposable ('nonDecomp = FALSE', default) or non-decomposable ('nonDecomp = TRUE')
    <li> The startup banner's model line now only distinguishes decomposable vs. non-decomposable graphs (the former "dense graph" banner line is gone along with 'fullCov')
    <li> Fixed a documentation mismatch in all four Scenario_*() simulators: the returned list element is named 'parSim', but every '@return'/'@examples' block referred to it as 'par' (this happened to still work via R's partial matching on '$', which is presumably why it went unnoticed). Docs and examples now correctly say 'parSim', matching the actual returned list element
    <li> Fixed a typo in Scenario_Dependence()'s '@return' doc where 'rho_X' was listed twice and 'rho_Y' was missing from the documented parameter list
    <li> Fixed Scenario_Mediation()'s '@return' summary tuple, which omitted the (correctly itemised and returned) 'mediation' element
    <li> Fixed MR2()'s '@return' doc, which omitted the (correctly returned) 'timeMR2' element
    <li> Fixed BMM_EM()'s '@param PIP' doc, which incorrectly said "returned by MR2" -- the mPIP is actually returned by 'PostProc()' (e.g. 'PostProc_output$gammaPost')
    <li> Added the missing '@return' item for BMM_EM()'s 'p' element (mixing proportions), which was listed in the summary tuple but never itemised
</ul>

Changes in **MR<sup>2</sup>** package version **0.1.2** (01/07/2025):
<ul>
    <li> Typos in the functions' descriptions have been corrected. Clarifications of some functions' outputs are also included
 </ul>   

Changes in **MR<sup>2</sup>** package version **0.1.1** (01/05/2025):
<ul>
    <li> All simulation functions now return '$betaHat_Y' and '$betaHat_X' to stress that summary-level data are generated
    <li> MR2() now considers explicitly as input 'betaHat_Y' and 'betaHat_Y' to stress that summary-level data are analysed
    <li> MR2() now also accepts a single value between (0, 1) for 'EVgamma'. It corresponds to the prior probability of exposure-outcome direct association. A single value is compulsory when only one exposure is employed
    <li> In MR2(), the option 'std'  has been removed. The analysis is automatically performed without the standardisation of the exposure summary-level data
    <li> Supporting routines for MR2() are now included in MR2() function
    <li> In BVS() (now included in the MR2() function), a bug in 'bPi <- hyperPar$bPi' assignment has been corrected
    <li> In PostProc(), rows and columns headers are added to all output matrices
</ul>
