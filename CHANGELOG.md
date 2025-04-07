
# News/Changelog 

Changes in **MR<sup>2</sup>** package version **1.1.0**:
<ul>
    <li> All simulation functions now return '$betaHat_Y' and '$betaHat_X' to stress that summary-level data are generated
    <li> MR2() now considers explicitely as input 'betaHat_Y' and 'betaHat_Y' to stress that summary-level data are analysed
    <li> MR2() now also accepts a single value between (0, 1) for 'EVgamma'. It corresponds to the prior probability of exposure-outcome causal association. A single value is compulsory when only one exposure is employed
    <li> In MR2(), the option 'std'  has been removed. The analysis is automatically performed without the standardization of the exposure summary-level data
    <li> Supporting routines for MR2() are now included in MR2() function
    <li> In BVS() (now included in the MR2() function), a bug in 'bPi <- hyperPar$bPi' assignment has been corrected
    <li> In PostProc(), rows and columns headers are added to all output matrices
</ul>


