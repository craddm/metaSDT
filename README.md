# metaSDT
A simple R package for calculating (meta-) SDT measures

This is an R implementation of Maniscalco and Lau's method of calculating metacognitive SDT measures using maximum likelihood estimation.

For further details on metacognitive/Type 2 SDT, see their website at http://www.columbia.edu/~bsm2105/type2sdt/

Currently, the fit_meta_d_MLE function is working - more commands are planned/partly implemented, but not yet available for use.

The fit_meta_d_MLE command requires data in the same format outlined on M & L's webpage, as follows:

Suppose there are two stimuli, A, and B, and three confidence ratings. The possible responses are

A1, A2, A3, B3, B2, B1.

Input to the function should be counts for each of these responses separately for each stimulus type.

fits <- fit_meta_d_MLE(nR_S1, nR_S2)

Output is a data frame with m-ratio etc.

devtools::install_github("craddm/metaSDT")

Note that the command is liable to change - although hopefully, the results won't!
