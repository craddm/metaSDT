# metaSDT
A simple R package for calculating (meta-) SDT measures

This is an R implementation of Maniscalco and Lau's method of calculating metacognitive SDT measures using maximum likelihood estimation.

For further details on metacognitive/Type 2 SDT, see their website at http://www.columbia.edu/~bsm2105/type2sdt/ and the following publications:
Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach for estimating metacognitive sensitivity from confidence ratings. Consciousness and Cognition, 21(1), 422–430. doi:10.1016/j.concog.2011.09.021

Maniscalco, B., & Lau, H. (2014). Signal detection theory analysis of type 1 and type 2 data: meta-d’, response-specific meta-d’, and the unequal variance SDT mode. In S. M. Fleming & C. D. Frith (Eds.), The Cognitive Neuroscience of Metacognition (pp.25-66). Springer.

Currently, the fit_meta_d_MLE function is working - more commands are planned/partly implemented, but not yet available for use.

If you use this function, cite the above papers and scripts on which it is based.

The fit_meta_d_MLE command requires data in the same format outlined on M & L's webpage, as follows:

Suppose there are two stimuli, A, and B, and three confidence ratings. The possible responses are

A1, A2, A3, B3, B2, B1.

Input to the function should be counts for each of these responses separately for each stimulus type.

fits <- fit_meta_d_MLE(nR_S1, nR_S2)

Output is a data frame with m-ratio etc.

devtools::install_github("craddm/metaSDT")

Note that the command is liable to change - although hopefully, the results won't!
