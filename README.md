
<!-- README.md is generated from README.Rmd. Please edit that file -->

# metaSDT

<!-- badges: start -->

[![R-CMD-check](https://github.com/craddm/metaSDT/workflows/R-CMD-check/badge.svg)](https://github.com/craddm/metaSDT/actions)
<!-- badges: end -->

This is an R implementation of Maniscalco and Lau’s methods of
calculating metacognitive SDT measures using maximum likelihood
estimation and minimization of the sum of squared errors.

## Installation

You can install `metaSDT`’s development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("craddm/metaSDT")
```

For further details on metacognitive/Type 2 SDT, see their website at
<http://www.columbia.edu/~bsm2105/type2sdt/> and the following
publications: Maniscalco, B., & Lau, H. (2012). A signal detection
theoretic approach for estimating metacognitive sensitivity from
confidence ratings. Consciousness and Cognition, 21(1), 422–430.
<doi:10.1016/j.concog.2011.09.021>

Maniscalco, B., & Lau, H. (2014). Signal detection theory analysis of
type 1 and type 2 data: meta-d’, response-specific meta-d’, and the
unequal variance SDT mode. In S. M. Fleming & C. D. Frith (Eds.), The
Cognitive Neuroscience of Metacognition (pp.25-66). Springer.

If you use these function, cite the above papers and scripts on which it
is based.

The fit\_meta\_d\_SSE and fit\_meta\_d\_MLE commands require data in the
same format outlined on M & L’s webpage, as follows:

Suppose there are two stimuli, A, and B, and three confidence ratings.
The possible responses are

A1, A2, A3, B3, B2, B1.

Input to the function should be counts for each of these responses
separately for each stimulus type.

So for example:

-   nR\_S1 &lt;- c(100,50,30,20,10,4)
-   nR\_S2 &lt;- c(4, 20, 21, 35, 60, 90)

fit\_MLE &lt;- fit\_meta\_d\_MLE(nR\_S1, nR\_S2) fit\_SSE &lt;-
fit\_meta\_d\_SSE(nR\_S1, nR\_S2)

Output is a data frame with m-ratio etc.
