#' Type 1 SDT for a 2AFC design
#'
#' This calculates standard type 1 SDT measures for 2AFC. The expected data frame format is one column indicating the stimulus (note that the first level of this factor will be treated as stimulus A), one column indicating the response (note - this should be coded with 1 as response A), and one column indicating the total number of responses of that type. Thus, there should be one row per combination of stimulus and response. If your data is in long format (i.e. one row per trial), you can use the \code{sdt_counts} function first to get the data into the expected format. The \code{type_1_sdt} function assumes that 1 = responded with first level of stimulus factor (e.g. 1 = stimulus A), and will calculate d-prime on that basis. Note that by default it adds a small constant to all cells to avoid boundary issues.
#'
#' @param df Data frame. See notes.
#' @param stimulus Column name for levels of the stimulus. Should be bare, unquoted. e.g. (stimulus = stimulus)
#' @param response Column name for responses. Should be bare, unquoted. e.g. (response = response).
#' @param counts Column name for totals. Should be bare, unquoted. Defaults to "total", as this is the column name output by \code{sdt_counts}.
#' @param s Ratio of standard deviations of stimulus types. Defaults to 1 (equal variance).
#' @param add_constant Adds a small constant to every cell to account for boundaries - i.e. log-linear correction. Default = TRUE.
#'
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#' @import dplyr
#' @import tidyr
#' @export

type_1_sdt <- function(df, stimulus = NULL, response = NULL, counts = total, s = 1, add_constant = TRUE) {
  stim_col <- enquo(stimulus)
  count_col <- enquo(counts)
  resp_col <- enquo(response)

  df <- group_by(df, !!stim_col)
  if (add_constant) {
    df <- mutate(df, proportions = ((!!count_col) + 1/n())/ (sum((!!count_col))+1))
  } else{
    df <- mutate(df, proportions = (!!count_col)/sum((!!count_col)))
  }

  reps_only <- filter(df, (!!resp_col) == 1)

  s1_HR <- reps_only$proportions[[1]]
  s1_FA <- reps_only$proportions[[2]]

  d_prime <- (1/s) * qnorm(s1_HR) - qnorm(s1_FA)
  c_raw <- (-1/(1+s)) * (qnorm(s1_HR)+qnorm(s1_FA))
  c_prime <- c_raw / d_prime
  data.frame(d_prime, c_raw, c_prime, s1_HR, s1_FA)
}

#' Convert trial-by-trial data to counts.
#'
#' This takes a trial-by-trial data frame and reduces it to counts. Intended mainly for use with the type 2 SDT fit_meta_d_MLE function. By default it will split the totals into multiple columns, one for each stimulus, with each row the total for a possible response. It is currently expected that confidence and response are combined into a single column  i.e. Response = "Definitely yes, maybe yes, maybe no, definitely no". Separate columns for confidence and response are not currently supported, but may be in the future.
#'
#' @param df Data frame containing single trial data.
#' @param stimulus Bare column name that contains the stimulus grouping of the trial (e.g. present versus absent).
#' @param response Bare column name that contains the response to be totalled over. (e.g. yes or no or any combination of confidence and response.)
#' @param split_resp Defaults to TRUE. Splits the counts into two columns, one for each stimulus.
#' @author Matt Craddock, \email{m.p.craddock@leeds.ac.uk}
#' @import dplyr
#' @import tidyr
#' @importFrom rlang UQE
#' @export

sdt_counts <- function(df, stimulus = NULL, response = NULL, split_resp = TRUE, ...) {
  stim_col <- enquo(stimulus)
  resp_col <- enquo(response)
  df <- group_by(df, !!stim_col, !!resp_col)
  df <- summarise(df, total = n())
  df <- ungroup(df)
  df <- complete_(df, c(rlang::UQE(stim_col), rlang::UQE(resp_col)), fill = list(total = 0))
  if (split_resp) {
    df <- spread_(df, quo_name(stim_col), "total")
  }
  return(df)
}

#' Fit Type 2 SDT using Maximum Likelihood Estimation.
#'
#' Provides a type-2 SDT analysis of data from a typical experiment in which observers discriminate between two response alternatives and provide ratings of confidence in their judgements.
#'
#' The expected input is two vectors, one for responses to each stimulus, encoding the observers response and confidence. For example, for two stimului labelled A and B, with three confidence ratings, participants could respond to stimulus A as follows:
#'
#' Response: A, rating: 3, count: 60
#'
#' Response: A, rating: 2, count: 30
#'
#' Response: A, rating: 1, count: 10
#'
#' Response: B, rating: 1, count: 7
#'
#' Response: B, rating: 2, count: 4
#'
#' Response: B, rating: 3, count: 1
#'
#' The appropriate vector would be nR_S1 <- c(60,30,10,7,4,1)
#'
#' For stimulus B, we would have the respective vector for responses to stimulus B, eg:
#'
#' Response: A, rating: 3, count: 4
#'
#' Response: A, rating: 2, count: 6
#'
#' Response: A, rating: 1, count: 11
#'
#' Response: B, rating: 1, count: 13
#'
#' Response: B, rating: 2, count: 23
#'
#' Response: B, rating: 3, count: 61
#'
#' nR_S2 <- c(4,6,11,13,23,61)
#'
#' The helper function \code{sdt_counts} can be used to get the data into the right format.
#'
#' The output is a data frame with various metacognitive measures, including m-ratio and meta-d, estimated using Maximum Likelihood Estimation. Currently, if more than 2 ratings are present in the data, the output will have multiple rows.
#'
#' For more details, see Maniscalco & Lau's webpage http://www.columbia.edu/~bsm2105/type2sdt/
#' Please cite that page and their articles if using this command.
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
#' @importFrom stats optim pnorm qnorm
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
  #Note that I suppress warnings from the optimizer about NaNs generated while calculating the log-likelihood. They don't effect the final estimates anyway, so the warnings are just annoying

  fit <- suppressWarnings(optim(par = guess,
               fit_meta_d_logL,
               gr = NULL,
               method = "BFGS",
               parameters = params,
               control = list("maxit" = 10000)))

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

  #Note that the S1units are not currently in the output - may change it to include them

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


#' Likelihood function for fitting meta-d.
#'
#' This is a likelihood function for use in MLE estimation, and shouldn't be called directly.
#'
#' @param x Starting guess for parameter values
#' @param parameters Various parameters such as the number of ratings, type 1 d-prime etc.
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

    prC_rS2[i] <- ((1 - pnorm(t2c1x[parameters$n_ratings + i], S2mu, S2sd)) - (1 - pnorm(t2c1x[parameters$n_ratings+ i+1], S2mu, S2sd))) / C_area_rS2
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


#' Function for calculating meta-d' by minimizing SSE
#'
#' Provides a type-2 SDT analysis of data from a typical experiment in which observers discriminate between two response alternatives and provide ratings of confidence in their judgements.
#'
#' Where fit_meta_d_MLE uses Maximum Likelihood Estimation, fit_meta_d_SSE works by finding the minimum sum of squared errors. As with the MLE method, input is expected as counts for each of two stimulus types.
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
#' The output is a dataframe with various metacognitive measures, including m-ratio and meta-d, estimated through minimizatoin of SSE.
#'
#'Currently, multiple rows will be returned when there are more than 2 confidence ratings.
#'
#' For more details, see Maniscalco & Lau's webpage http://www.columbia.edu/~bsm2105/type2sdt/
#' Please cite that page and their articles if using this command.
#'
#'
#'@import dplyr
#'@import purrr
#'@export

fit_meta_d_SSE <- function(nR_S1, nR_S2, s = 1, d_min = -5, d_max = 5, d_grain = .01, add_constant = TRUE) {

  if (add_constant) {
    nR_S1 <- nR_S1 + (1/length(nR_S1))
    nR_S2 <- nR_S2 + (1/length(nR_S2))
  }

  n_ratings <- length(nR_S1) / 2
  n_criteria <- 2 * n_ratings - 1

  S1_HR <- sum(nR_S2[(n_ratings+1):(n_ratings*2)])/sum(nR_S2)
  S1_FA <- sum(nR_S1[(n_ratings+1):(n_ratings*2)])/sum(nR_S1)

  d_1 <- (1/s) * qnorm(S1_HR) - qnorm(S1_FA)
  c_1 <- (-1 / (1 + s)) * (qnorm(S1_HR) + qnorm(S1_FA))
  c_prime <- c_1 / d_1

  obs_HR2_rS1 <- NULL
  obs_FAR2_rS1 <- NULL

  for (i in 1:(n_ratings-1)) {
    obs_HR2_rS1[i] <- sum(nR_S1[1:i]) / sum(nR_S1[1:n_ratings])
    obs_FAR2_rS1[i] <- sum(nR_S2[1:i]) / sum(nR_S2[1:n_ratings])
  }

  obs_HR2_rS2 <- NULL
  obs_FAR2_rS2 <- NULL
  for (i in (n_ratings+2):(2*n_ratings)) {
    obs_HR2_rS2[(i - (n_ratings+2) + 1)] <- sum(nR_S2[i:(n_ratings*2)]) / sum(nR_S2[(n_ratings+1):(n_ratings*2)])
    obs_FAR2_rS2[(i - (n_ratings+2) + 1)] <- sum(nR_S1[i:(n_ratings*2)]) / sum(nR_S1[(n_ratings+1):(n_ratings*2)])
  }

  d_grid <- seq(d_min, d_max, by = d_grain)
  c_grid <- c_prime * d_grid

  S1mu <- -d_grid / 2
  S2mu <- d_grid / 2
  S1sd <- 1
  S2sd <- 1 / s

  bounds <- 5 * max(S1sd, S2sd)
  SSEmin <- Inf

  param_space <- map2(S1mu, S2mu, ~seq(.x - bounds, .y + bounds, by = .001)) # param space

  min_z <- map2(param_space, c_grid, ~min(abs(.x - .y)))
  c_ind <- map2(param_space, c_grid, ~which.min(abs(.x - .y)))

  HRs <- map2(param_space, S2mu,
              ~ 1 - (pnorm(.x, .y, S2sd)))
  FARs <- map2(param_space, S1mu,
               ~ 1 - (pnorm(.x, .y, S1sd)))

  # fit type 2 data for S1 responses
  est_HR2s_rS1 <- map2(FARs, c_ind, ~(1 - .x[1:.y]) / (1 - .x[.y]))
  est_FAR2s_rS1 <- map2(HRs, c_ind, ~(1 - .x[1:.y]) / (1 - .x[.y]))

  SSE <- NULL
  SSE_rS1 <- matrix(data = NA, nrow = length(d_grid), ncol = n_ratings-1)
  rS1_ind <- matrix(data = NA, nrow = length(d_grid), ncol = n_ratings-1)
  for (n in (1:(n_ratings-1))) {
    SSE <- map2(est_HR2s_rS1, est_FAR2s_rS1,
              ~(.x - obs_HR2_rS1[[n]]) ^2 + (.y - obs_FAR2_rS1[[n]]) ^2)
    SSE_rS1[ ,n] <- map_dbl(SSE, min)
    rS1_ind[ ,n] <- map_dbl(SSE, which.min)
  }

  # fit type 2 data for S2 responses
  est_HR2s_rS2 <- map2(HRs, c_ind, ~ .x[.y:length(.x)] / .x[.y])
  est_FAR2s_rS2 <- map2(FARs, c_ind, ~ .x[.y:length(.x)] / .x[.y])

  SSE_rS2 <- matrix(data = NA, nrow = length(d_grid), ncol = n_ratings-1)
  rS2_ind <- matrix(data = NA, nrow = length(d_grid), ncol = n_ratings-1)
  for (n in (1:(n_ratings-1))) {
    SSE <- map2(est_HR2s_rS2, est_FAR2s_rS2,
                ~(.x - obs_HR2_rS2[[n]]) ^2 + (.y - obs_FAR2_rS2[[n]]) ^2)
    SSE_rS2[ ,n] <- map_dbl(SSE, min)
    rS2_ind[ ,n] <- map_dbl(SSE, which.min)
  }

  # update analysis
  SSEtot <- rowSums(SSE_rS1) + rowSums(SSE_rS2)
  SSEmin <- Inf

  meta_d <- vector("numeric", 1)
  meta_c <- vector("numeric", 1)
  t2c_rS1 <- NULL
  t2c_rS2 <- NULL

  min_SSE <- which.min(SSEtot)
  meta_d <- d_grid[[min_SSE]]
  meta_c <- c_grid[[min_SSE]]
  t2c_rS1 <- param_space[[min_SSE]][rS1_ind[[min_SSE]]]
  t2c_rS2 <- param_space[[min_SSE]][c_ind[[min_SSE]] + rS2_ind[[min_SSE]] - 1]
  est_HR2_rS1  <- est_HR2s_rS1[[min_SSE]][rS1_ind[[min_SSE]]]
  est_FAR2_rS1 <- est_FAR2s_rS1[[min_SSE]][rS1_ind[[min_SSE]]]
  est_HR2_rS2  <- est_HR2s_rS2[[min_SSE]][rS2_ind[[min_SSE]]]
  est_FAR2_rS2 <- est_FAR2s_rS2[[min_SSE]][rS2_ind[[min_SSE]]]

  out <- data.frame(da = d_1,
                    meta_da = meta_d,
                    M_diff = meta_d - d_1,
                    M_ratio = meta_d / d_1,
                    meta_ca = meta_c,
                    s = s,
                    t2ca_rS1 = t2c_rS1,
                    t2ca_rS2 = t2c_rS2,
                    SSE = SSEtot[[min_SSE]],
                    est_HR2_rS1 = est_HR2_rS1,
                    obs_HR2_rS1 = obs_HR2_rS1,
                    est_HR2_rS2 = est_HR2_rS2,
                    obs_HR2_rS2 = obs_HR2_rS2,
                    est_FAR2_rS1 = est_FAR2_rS1,
                    obs_FAR2_rS1 = obs_FAR2_rS1,
                    est_FAR2_rS2 = est_FAR2_rS2,
                    obs_FAR2_rS2 = obs_FAR2_rS2)
  return(out)
}
