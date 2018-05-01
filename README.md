# metaSDT

A simple R package for calculating (meta-) SDT measures

The package implements several different methods of calculating Type 1 and Type 2 Signal Detection measures.

To install, use the following command:

`devtools::install_github("craddm/metaSDT")`

## Type 1 SDT 

Two functions exist for calculating Type 1 SDT measures from count data.

- `type_1_sdt()`
    - For use in 2AFC tasks (i.e. response is Stimulus A or Stimulus B)
- `rating_sdt()`
    - For use when responses are on a rating scale (e.g. 1-6, with 1 being Definitely Stimulus A and 6 being Definitely Stimulus B)

## Type 2 SDT

The package implements several different measures of metacognitive sensitivity / Type 2 SDT.

Maniscalco and Lau's methods of calculating metacognitive SDT measures using maximum likelihood estimation and minimization of the sum of squared errors are implemented in two separate functions.

- `fit_meta_d_MLE()`
- `fit_meta_d_sse()`

For further details on these methods, see their website at http://www.columbia.edu/~bsm2105/type2sdt/ and the following publications:

Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach for estimating metacognitive sensitivity from confidence ratings. Consciousness and Cognition, 21(1), 422–430. doi:10.1016/j.concog.2011.09.021

Maniscalco, B., & Lau, H. (2014). Signal detection theory analysis of type 1 and type 2 data: meta-d’, response-specific meta-d’, and the unequal variance SDT mode. In S. M. Fleming & C. D. Frith (Eds.), The Cognitive Neuroscience of Metacognition (pp.25-66). Springer.

If you use these functions, cite the above papers and scripts on which it is based.

A third method of calculating metacognitive measures is metad' *balance* (Barrett, Dienes, & Seth, 2013)

- `fit_meta_d_bal()`

See:
Barrett, A.B., Dienes, Z., & Seth, A.K. (2013). Measures of metacognition on signal-detection theoretic models. Psychological Methods, 18(4), 535-552. doi: 10.1037/a0033268

### Usage guidelines

The `fit_meta_d_SSE`, `fit_meta_d_MLE`, and `fit_meta_d_bal` commands require data in the format outlined on M & L's webpage, which is as follows:

Suppose there are two stimuli, A, and B, and three confidence ratings. The possible responses are

1. A1
2. A2
3. A3
4. B3
5. B2
6. B1

Input to the function should be counts for each of these responses separately for each stimulus type.

So for example:

- nR_S1 <- c(100, 50, 30, 20, 10, 4)
- nR_S2 <- c(4, 20, 21, 35, 60, 90)

`fit_MLE <- fit_meta_d_MLE(nR_S1, nR_S2)`
`fit_SSE <- fit_meta_d_SSE(nR_S1, nR_S2)`
`fit_bal <- fit_meta_d_bal(nR_S1, nR_S2)`

Output is a data frame with m-ratio etc.

## Generating count data

If you have a dataframe with one row per observation/trial, you can use the `sdt_counts()` function to generate a summary of counts for each combination of trial type (e.g. stimulus present vs absent) and response (e.g. present or absent). Note that by default this will split the counts into two columns, one for each trial type, with each row being a response.
