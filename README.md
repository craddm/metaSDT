
<!-- README.md is generated from README.Rmd. Please edit that file -->

# metaSDT

<!-- badges: start -->

[![R-CMD-check](https://github.com/craddm/metaSDT/workflows/R-CMD-check/badge.svg)](https://github.com/craddm/metaSDT/actions)
[![DOI](https://zenodo.org/badge/99712128.svg)](https://zenodo.org/badge/latestdoi/99712128)
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

For further details on metacognitive/Type 2 SDT, see Maniscalco & Lau's website at
<http://www.columbia.edu/~bsm2105/type2sdt/> and the following
publications:

Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach
for estimating metacognitive sensitivity from confidence ratings.
Consciousness and Cognition, 21(1), 422–430.
<doi:10.1016/j.concog.2011.09.021>

Maniscalco, B., & Lau, H. (2014). Signal detection theory analysis of
type 1 and type 2 data: meta-d’, response-specific meta-d’, and the
unequal variance SDT mode. In S. M. Fleming & C. D. Frith (Eds.), The
Cognitive Neuroscience of Metacognition (pp.25-66). Springer.

If you use these functions, cite the above papers and scripts on which
it is based.

The `fit_meta_d_SSE` and `fit_meta_d_MLE` commands require data in the
same format outlined on M & L’s webpage, as follows:

Suppose there are two stimuli, A, and B, and three confidence ratings.
The possible responses are

A1, A2, A3, B3, B2, B1.

Input to the function should be counts for each of these responses
separately for each stimulus type.

So for example:

``` r
library(metaSDT)
nR_S1 <- c(100, 50, 30, 20, 10, 4)
nR_S2 <- c(4, 20, 21, 35, 60, 90)
fit_MLE <- fit_meta_d_MLE(nR_S1,
                          nR_S2)
fit_MLE
#>         da s meta_da     M_diff  M_ratio    meta_ca   t2ca_rS1  t2ca_rS2
#> 1 1.845043 1 1.87284 0.02779706 1.015066 0.07108153 -1.0297441 0.5524547
#> 2 1.845043 1 1.87284 0.02779706 1.015066 0.07108153 -0.3710477 1.2052824
#>        logL est_HR2_rS1 est_HR2_rS2 est_FAR2_rS1 est_FAR2_rS2 obs_HR2_rS1
#> 1 -444.1742   0.8469277   0.8052557    0.4938669    0.4352026   0.8328717
#> 2 -444.1742   0.5489194   0.4885083    0.1273860    0.1027034   0.5549400
#>   obs_HR2_rS2 obs_FAR2_rS1 obs_FAR2_rS2       t1c1
#> 1   0.8104223   0.53479853    0.4154589 0.07002653
#> 2   0.4860737   0.09157509    0.1207729 0.07002653
fit_SSE <- fit_meta_d_SSE(nR_S1,
                          nR_S2)
fit_SSE
#>         da meta_da     M_diff  M_ratio    meta_ca s t2ca_rS1 t2ca_rS2
#> 1 1.845043    1.89 0.04495688 1.024366 0.07173282 1   -1.057    0.568
#> 2 1.845043    1.89 0.04495688 1.024366 0.07173282 1   -1.057    0.568
#>          SSE est_HR2_rS1 obs_HR2_rS1 est_HR2_rS2 obs_HR2_rS2 est_FAR2_rS1
#> 1 0.00246936   0.5386789   0.5549400   0.7999733   0.8104223     0.118341
#> 2 0.00246936   0.5386789   0.8328717   0.7999733   0.4860737     0.118341
#>   obs_FAR2_rS1 est_FAR2_rS2 obs_FAR2_rS2       t1c1
#> 1   0.09157509    0.4214078    0.4154589 0.07002653
#> 2   0.53479853    0.4214078    0.1207729 0.07002653
```

Output is a data frame with m-ratio etc.
