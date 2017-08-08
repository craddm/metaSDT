#' Type 1 SDT for a 2AFC design
#'
#' This calculates standard, type 1 SDT measures.
#'
#' @param df Raw data frame in long format. Assumes counts already calculated. Expects the response column to be 1 or 0 for reported yes or reported no.
#' @param stimulus Column name for levels of the stimulus. Should be bare, unquoted. e.g. (stimulus = stimulus)
#' @param response Col
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#' @import dplyr
#' @import tidyr

type_1_sdt <- function(df, stimulus = NULL, response = NULL, counts = total, confidence = NULL, s = 1, add_constant = TRUE) {
  stim_col <- enquo(stimulus)
  count_col <- enquo(counts)
  resp_col <- enquo(response)

  df <- group_by(df, !!stim_col)
  if (add_constant) {
    df <- mutate(df, proportions = ((!!count_col) + 1/n())/ (sum((!!count_col))+1))
  } else{
    df <- mutate(df, proportions = (!!count_col)/sum((!!count_col)))
  }

  s1_HR <- filter(df, (!!resp_col) == 1)

  s1_FA <- sum(df$proportions[5:6])
  d_prime <- (1/s) * qnorm(s1_HR) - qnorm(s1_FA)
  c_raw <- (-1/(1+s)) * (qnorm(s1_HR)+qnorm(s1_FA))
  c_prime <- c_raw / d_prime
  data.frame(d_prime, c_raw, c_prime, s1_HR, s1_FA)
}

#' Convert to counts
#'
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#' @import dplyr
#' @import tidyr

sdt_counts <- function(df, stimulus = NULL, response = NULL, ...) {
  stim_col <- enquo(stimulus)
  resp_col <- enquo(response)
  df <- group_by(df, !!stim_col, !!resp_col)
  df <- summarise(df, total = n())
  df <- ungroup(df)
  df <- complete_(df, c(rlang::UQE(stim_col), rlang::UQE(resp_col)), fill = list(total = 0))
}

#' Fit Type 2 SDT using Maximum Likelihood Estimation.
#'
#' Provides a type-2 SDT analysis of data from a typical experiment in which observers discriminate between two response alternatives and provide ratings of confidence in their judgements.
#'
#' The expected input is two vectors, one for responses to each stimulus, encoding the observers response and confidence. For example, for two stimului labelled A and B, with three confidence ratings, participants could respond to stimulus A as follows:
#' Response: A, rating: 3, count: 60
#' Response: A, rating: 2, count: 30
#' Response: A, rating: 1, count: 10
#' Response: B, rating: 1, count: 7
#' Response: B, rating: 2, count: 4
#' Response: B, rating: 3, count: 1
#'
#' The appropriate vector would be nR_S1 <- c(60,30,10,7,4,1)
#'
#' For stimulus B, we would have the respective vector for responses to stimulus B, eg:
#' Response: A, rating: 3, count: 4
#' Response: A, rating: 2, count: 6
#' Response: A, rating: 1, count: 11
#' Response: B, rating: 1, count: 13
#' Response: B, rating: 2, count: 23
#' Response: B, rating: 3, count: 61
#'
#' nR_S2 <- c(4,6,11,13,23,61)
#'
#' The output is a dataframe with various metacognitive measures, including m-ratio and meta-d, estimated using Maximum Likelihood Estimation.
#'
#' For more details, see Maniscalco & Lau's webpage http://www.columbia.edu/~bsm2105/type2sdt/
#'
#' @param nR_S1 Responses to S1 stimulus. See below for advice.
#' @param nR_S2 Responses to S2 stimulus. See below for advice.
#' @param s Ratio of standard deviations for the S1 and S2 stimulus. Defaults to 1.
#' @param add_constant Adds a small constant to the data (1/number of possible responses) to account for 0 or 1 values. Defaults to TRUE for ease of use across multiple datasets.
#'
#' @author Maniscalco & Lau. Ported to R by Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#' @references Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach for estimating metacognitive sensitivity from confidence ratings. Consciousness and Cognition. http://dx.doi.org/10.1016/j.concog.2011.09.021
#' @import dplyr
#' @import tidyr
#' @export
#'

fit_meta_d_MLE <- function(nR_S1, nR_S2, s = 1, add_constant = TRUE) {

  if (add_constant) {
    nR_S1 <- nR_S1 + (1/length(nR_S1))
    nR_S2 <- nR_S2 + (1/length(nR_S2))
  }

  n_ratings <- length(nR_S1) / 2
  n_criteria <- 2 * n_ratings - 1

  A <- matrix(0,nrow = n_criteria - 2, ncol = n_criteria)
  for (i in 2:n_criteria - 1) {
    A[i-1, i:(i + 1)] <- c(1, -1)
  }

  b <- rep(-1e-5, n_criteria - 2)

  #LB <- c(-10, -20*(rep(1, (n_criteria-1)/2)), rep(0, (n_criteria-1)/2))
  #UB <- c(10, rep(0, (n_criteria-1)/2), 20*(rep(1, (n_criteria-1)/2)))

  ## set up initial guess at parameter values
  ratingHR  <- NULL
  ratingFAR <- NULL
  for (c in 2:(n_ratings * 2)) {
    ratingHR[c-1] <- sum(nR_S2[c:length(nR_S2)]) / sum(nR_S2)
    ratingFAR[c-1] <- sum(nR_S1[c:length(nR_S1)]) / sum(nR_S1)
  }

  t1_index <- n_ratings
  t2_index <- setdiff(1:(2 * n_ratings - 1), t1_index)

  d1 <- (1 / s) * qnorm(ratingHR[t1_index]) - qnorm(ratingFAR[t1_index])
  meta_d1 <- d1

  c1 <- (-1/(1 + s)) * ( qnorm(ratingHR) + qnorm(ratingFAR) )
  t1c1 <- c1[t1_index]
  t2c1 <- c1[t2_index]

  guess <- c(meta_d1, t2c1 - (meta_d1 * t1c1 / d1))

  params <- list("n_ratings" = n_ratings,
                 "d1" = d1,
                 "nR_S1" = nR_S1,
                 "nR_S2" = nR_S2,
                 "t1c1" = t1c1,
                 "s" = s)

  fit <- suppressWarnings(optim(par = guess,
               fit_meta_d_logL,
               gr = NULL,
               method = "BFGS",
               parameters = params))

  meta_d1 <- fit$par[[1]]
  t2c1 <- fit$par[2:length(fit$par)] + (meta_d1 * t1c1 / d1)
  logL <- -fit$value

  I_nR_rS2 <- nR_S1[(n_ratings + 1):length(nR_S1)]
  I_nR_rS1 <- nR_S2[seq(n_ratings, 1, -1)]

  C_nR_rS2 <- nR_S2[(n_ratings+1):length(nR_S2)]
  C_nR_rS1 <- nR_S1[seq(n_ratings, 1, -1)]
  obs_FAR2_rS2 <- matrix(0, nrow = n_ratings - 1)
  obs_HR2_rS2 <- matrix(0, nrow = n_ratings - 1)
  obs_FAR2_rS1 <- matrix(0, nrow = n_ratings - 1)
  obs_HR2_rS1 <- matrix(0, nrow = n_ratings - 1)

  for (i in 2:n_ratings) {
    obs_FAR2_rS2[i - 1] <- sum(I_nR_rS2[i:length(I_nR_rS2)] ) / sum(I_nR_rS2)
    obs_HR2_rS2[i - 1]  <- sum(C_nR_rS2[i:length(I_nR_rS2)] ) / sum(C_nR_rS2)

    obs_FAR2_rS1[i-1] <- sum(I_nR_rS1[i:length(I_nR_rS2)] ) / sum(I_nR_rS1)
    obs_HR2_rS1[i-1]  <- sum(C_nR_rS1[i:length(I_nR_rS2)] ) / sum(C_nR_rS1)
  }

  S1mu <- -meta_d1 / 2
  S1sd <- 1
  S2mu <-  meta_d1 / 2
  S2sd <- S1sd / s
  mt1c1 <- (meta_d1 * t1c1 / d1);

  C_area_rS2 <- 1 - pnorm(mt1c1,S2mu,S2sd);
  I_area_rS2 <- 1 - pnorm(mt1c1,S1mu,S1sd);
  C_area_rS1 <- pnorm(mt1c1,S1mu,S1sd);
  I_area_rS1 <- pnorm(mt1c1,S2mu,S2sd);

  est_FAR2_rS2 <- NULL
  est_FAR2_rS1 <- NULL
  est_HR2_rS2 <- NULL
  est_HR2_rS1 <- NULL

  for (i in 1:(n_ratings - 1)) {

    t2c1_lower <- t2c1[n_ratings - i]
    t2c1_upper <- t2c1[n_ratings - 1 + i]

    I_FAR_area_rS2 <- 1-pnorm(t2c1_upper, S1mu, S1sd)
    C_HR_area_rS2  <- 1-pnorm(t2c1_upper, S2mu, S2sd)

    I_FAR_area_rS1 <- pnorm(t2c1_lower, S2mu, S2sd)
    C_HR_area_rS1  <- pnorm(t2c1_lower, S1mu, S1sd)


    est_FAR2_rS2[i] <- I_FAR_area_rS2 / I_area_rS2
    est_HR2_rS2[i]  <- C_HR_area_rS2 / C_area_rS2

    est_FAR2_rS1[i] <- I_FAR_area_rS1 / I_area_rS1
    est_HR2_rS1[i]  <- C_HR_area_rS1 / C_area_rS1
  }

  da <- sqrt(2 / (1 + s^2)) * s * d1
  meta_da <- sqrt(2 / (1 + s^2)) * s * meta_d1
  mt1c1 <- (meta_d1 * t1c1 / d1)
  t2ca <- (sqrt(2) * s / sqrt(1 + s ^2)) * t2c1
  S1units <- data.frame(d1 = d1,
                        meta_d1 = meta_d1,
                        s = s,
                        meta_c1 = mt1c1,
                        t2c1_rS1 = t2c1[1:n_ratings-1],
                        t2c1_rS2 = t2c1[n_ratings:length(t2c1)])

  fit <- data.frame(da = da,
                    s = s,
                    meta_da = meta_da,
                    M_diff = meta_da - da,
                    M_ratio = meta_da / da,
                    meta_ca = (sqrt(2) * s / sqrt(1 + s ^2) * mt1c1),
                    t2ca_rS1 = t2ca[1:n_ratings-1],
                    t2ca_rS2 = t2ca[n_ratings:length(t2ca)],
                    #          S1units = list(S1units),
                    logL = logL,
                    est_HR2_rS1 = est_HR2_rS1,
                    est_HR2_rS2 = est_HR2_rS2,
                    est_FAR2_rS1 = est_FAR2_rS1,
                    est_FAR2_rS2 = est_FAR2_rS2,
                    obs_HR2_rS1 = obs_HR2_rS1,
                    obs_HR2_rS2 = obs_HR2_rS2,
                    obs_FAR2_rS1 = obs_FAR2_rS1,
                    obs_FAR2_rS2 = obs_FAR2_rS2

  )
  return(fit)
}


#' Likelihood function for fitting meta-d
#'
#' @author Maniscalco and Lau. Ported to R by Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#'

fit_meta_d_logL <- function(x, parameters) {
  meta_d1 <- x[1]
  t2c1 <- x[2:length(x)]

  S1mu <- -meta_d1/2
  S1sd <- 1
  S2mu <- meta_d1/2
  S2sd <- S1sd/parameters$s

  S1mu <- S1mu - (meta_d1 * (parameters$t1c1 / parameters$d1))
  S2mu <- S2mu - (meta_d1 * (parameters$t1c1 / parameters$d1))

  t1c1 <- 0
  nC_rS1 <- matrix(0, ncol = parameters$n_ratings)
  nI_rS1 <- matrix(0, ncol = parameters$n_ratings)
  nC_rS2 <- matrix(0, ncol = parameters$n_ratings)
  nI_rS2 <- matrix(0, ncol = parameters$n_ratings)

  for (i in 1:parameters$n_ratings) {

    # S1 responses
    nC_rS1[i] <- parameters$nR_S1[i]
    nI_rS1[i] <- parameters$nR_S2[i]

    # S2 responses
    nC_rS2[i] <- parameters$nR_S2[parameters$n_ratings + i]
    nI_rS2[i] <- parameters$nR_S1[parameters$n_ratings + i]
  }

  # get type 2 probabilities
  C_area_rS1 <- pnorm(t1c1,S1mu,S1sd)
  I_area_rS1 <- pnorm(t1c1,S2mu,S2sd)

  C_area_rS2 <- 1 - pnorm(t1c1,S2mu,S2sd)
  I_area_rS2 <- 1 - pnorm(t1c1,S1mu,S1sd)

  t2c1x <- c(-Inf, t2c1[1:parameters$n_ratings - 1], t1c1, t2c1[parameters$n_ratings:length(t2c1)], Inf)
  prC_rS1 <- matrix(0, ncol = parameters$n_ratings)
  prI_rS1 <- matrix(0, ncol = parameters$n_ratings)
  prC_rS2 <- matrix(0, ncol = parameters$n_ratings)
  prI_rS2 <- matrix(0, ncol = parameters$n_ratings)

  for (i in 1:parameters$n_ratings) {
    prC_rS1[i] <- (pnorm(t2c1x[i + 1], S1mu, S1sd) - pnorm(t2c1x[i], S1mu, S1sd) ) / C_area_rS1
    prI_rS1[i] <- (pnorm(t2c1x[i + 1], S2mu, S2sd) - pnorm(t2c1x[i], S2mu, S2sd) ) / I_area_rS1

    prC_rS2[i] <- ((1 - pnorm(t2c1x[parameters$n_ratings + i], S2mu, S2sd)) - (1-pnorm(t2c1x[parameters$n_ratings+ i+1], S2mu, S2sd))) / C_area_rS2
    prI_rS2[i] <- ((1 - pnorm(t2c1x[parameters$n_ratings + i], S1mu, S1sd)) - (1 - pnorm(t2c1x[parameters$n_ratings+i+1], S1mu, S1sd))) / I_area_rS2
  }


  # calculate log likelihood
  logL <- 0
  for (i in 1:parameters$n_ratings) {
    logL <- logL + nC_rS1[i] * log(prC_rS1[i]) + nI_rS1[i] * log(prI_rS1[i]) +
      nC_rS2[i] * log(prC_rS2[i]) + nI_rS2[i] * log(prI_rS2[i])
  }

  if (is.nan(logL)) {
    logL <- -Inf
  }
  logL <- -logL
  return(logL)
}
