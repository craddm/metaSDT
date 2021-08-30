#' Type 1 SDT for a 2AFC design
#'
#' Calculate standard type 1 SDT measures for 2AFC.
#'
#' The expected data frame format is one column indicating the stimulus (note
#' that the first level of this factor will be treated as stimulus A), one
#' column indicating the response (note - this should be coded with 1 as
#' response A), and one column indicating the total number of responses of that
#' type. Thus, there should be one row per combination of stimulus and response.
#' If your data is in long format (i.e. one row per trial), you can use the
#' \code{sdt_counts} function first to get the data into the expected format.
#' The \code{type_1_sdt} function assumes that 1 = responded with first level of
#' stimulus factor (e.g. 1 = stimulus A), and will calculate d-prime on that
#' basis. Note that by default it adds a small constant to all cells to avoid
#' boundary issues.
#'
#' @param df Data frame. See notes.
#' @param stimulus Column name for levels of the stimulus. Should be bare,
#'   unquoted. Default is `stimulus`
#' @param response Column name for responses. Should be bare, unquoted. Default
#'   is `response`.
#' @param counts Column name for totals. Should be bare, unquoted. Defaults to
#'   "total", as this is the column name output by \code{sdt_counts}.
#' @param s Ratio of standard deviations of stimulus types. Defaults to 1 (equal
#'   variance).
#' @param add_constant Adds a small constant to every cell to account for
#'   boundaries - i.e. log-linear correction. Default = TRUE.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#' @import dplyr
#' @import tidyr
#' @family type_1_sdt sdt_counts
#' @examples
#' type_1_test <- data.frame(
#'                  expand.grid(
#'                    stimulus = c("A", "B"),
#'                    response = c(1, 2)
#'                  ),
#'                  total = c(30, 10, 8, 32)
#'                  )
#' type_1_sdt(type_1_test)
#' @export

type_1_sdt <- function(df,
                       stimulus = stimulus,
                       response = response,
                       counts = total,
                       s = 1,
                       add_constant = TRUE) {

  stim_col <- dplyr::enquo(stimulus)
  count_col <- dplyr::enquo(counts)
  resp_col <- dplyr::enquo(response)

  df <- dplyr::group_by(df, !!stim_col)

  if (add_constant) {
    df <- dplyr::mutate(df, proportions = ( (!!count_col) + 1 /
                                       n()) / (sum( (!!count_col)) + 1))
  } else{
    df <- dplyr::mutate(df, proportions = (!!count_col) / sum( (!!count_col)))
  }

  reps_only <- dplyr::filter(df, (!!resp_col) == 1)

  s1_HR <- reps_only$proportions[[1]]
  s1_FA <- reps_only$proportions[[2]]

  d_prime <- (1 / s) * qnorm(s1_HR) - qnorm(s1_FA)
  c_raw <- (-1 / (1 + s)) * (qnorm(s1_HR) + qnorm(s1_FA))
  c_prime <- c_raw / d_prime
  data.frame(d_prime,
             c_raw,
             c_prime,
             s1_HR,
             s1_FA)
}

#' Type 1 SDT for rating data
#'
#' Provides a Type-1 SDT analysis of data from a rating experiment in which
#' observers discriminate between two response alternatives and provide ratings
#' of confidence in their judgements.
#'
#' The expected input is two vectors, one for responses to each stimulus,
#' encoding the observers response and confidence. For example, for two stimului
#' labelled A and B, with three confidence ratings, participants could respond
#' to stimulus A as follows:
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
#' For stimulus B, we would have the respective vector for responses to stimulus
#' B, eg:
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
#' The helper function \code{sdt_counts} can be used to get the data into the
#' right format.
#'
#' The output is a data frame with various Type-1 measures, one for each
#' criterion.
#'
#' @param nr_s1 Responses to S1 stimulus. See below for advice.
#' @param nr_s2 Responses to S2 stimulus. See below for advice.
#' @param add_constant Adds a small constant to every cell to account for
#'   boundaries - i.e. log-linear correction. Default = TRUE.
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @examples
#' nr_s1 <- c(60,30,10,7,4,1)
#' nr_s2 <- c(4,6,11,13,23,61)
#' rating_sdt(nr_s1, nr_s2)
#' @export

rating_sdt <- function(nr_s1,
                       nr_s2,
                       add_constant = FALSE) {

  if (add_constant) {
    nr_s1 <- nr_s1 + (1 / length(nr_s1))
    nr_s2 <- nr_s2 + (1 / length(nr_s2))
  }

  n_ratings <- length(nr_s1)
  n_criteria <- n_ratings - 1

  false_alarms <- NULL
  hit_rates <- NULL

  for (i in 1:(n_criteria)) {
    hit_rates[i] <- sum(nr_s1[1:i]) / sum(nr_s1)
    false_alarms[i] <- sum(nr_s2[1:i]) / sum(nr_s2)
  }

  z_fa <- qnorm(false_alarms)
  z_hr <- qnorm(hit_rates)
  dprime <- z_hr - z_fa
  Adprime <- pnorm(dprime / sqrt(2))

  diff_hr <- hit_rates - false_alarms

  Aprime <- 0.5 + (abs(diff_hr) / (diff_hr)) * ((diff_hr) ^ 2 + abs(diff_hr)) /
    (4 * max(hit_rates, false_alarms) - 4 * hit_rates * false_alarms)

  data.frame(false_alarms, hit_rates, z_fa, z_hr, dprime, Adprime, Aprime)
}


#' Convert trial-by-trial data to counts.
#'
#' Reduces a trial-by-trial data frame to counts.
#'
#' Intended mainly for use with the type 2 SDT fit_meta_d_MLE, fit_meta_d_SSE,
#' and fit_meta_d_bal functions. By default it will split the totals into
#' multiple columns, one for each stimulus, with each row the total for a
#' possible response. It is currently expected that confidence and response are
#' combined into a single column  i.e. Response = "Definitely yes, maybe yes,
#' maybe no, definitely no". Separate columns for confidence and response are
#' not currently supported, but may be in the future.
#'
#' @param df Data frame containing single trial data.
#' @param stimulus Bare column name that contains the stimulus grouping of the
#'   trial (e.g. present versus absent).
#' @param response Bare column name that contains the response to be totalled
#'   over. (e.g. yes or no or any combination of confidence and response.)
#' @param split_resp Defaults to TRUE. Splits the counts into two columns, one
#'   for each stimulus.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#' @import dplyr
#' @import tidyr
#' @import rlang
#' @family type_1_sdt fit_meta_d_MLE
#' @export

sdt_counts <- function(df, stimulus = NULL, response = NULL,
                       split_resp = TRUE) {
  stim_col <- rlang::enquo(stimulus)
  resp_col <- rlang::enquo(response)
  df <- dplyr::group_by(df, !!stim_col, !!resp_col)
  df <- dplyr::summarise(df, total = n())
  df <- dplyr::ungroup(df)
  df <- tidyr::complete(df, !!stim_col, !!resp_col, fill = list(total = 0))
  if (split_resp) {
    df <- tidyr::spread(df, !!stim_col, "total")
  }
  df
}

#' Fit Type 2 SDT using Maximum Likelihood Estimation.
#'
#' Provides a type-2 SDT analysis of data from a typical experiment in which
#' observers discriminate between two response alternatives and provide ratings
#' of confidence in their judgements.
#'
#' The expected input is two vectors, one for responses to each stimulus,
#' encoding the observers response and confidence. For example, for two stimului
#' labelled A and B, with three confidence ratings, participants could respond
#' to stimulus A as follows:
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
#' For stimulus B, we would have the respective vector for responses to stimulus
#' B, eg:
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
#' The helper function \code{sdt_counts} can be used to get the data into the
#' right format.
#'
#' The output is a data frame with various metacognitive measures, including
#' m-ratio and meta-d, estimated using Maximum Likelihood Estimation. The output
#' will have n-1 rows, where n is the number of confidence ratings. Most
#' measures are duplicated across rows, but some vary, as they reflect cutpoints
#' between the confidence levels.
#'
#' For more details, see Maniscalco & Lau's webpage
#' http://www.columbia.edu/~bsm2105/type2sdt/ Please cite that page and their
#' articles if using this command.
#'
#' @param nR_S1 Responses to S1 stimulus. See below for advice.
#' @param nR_S2 Responses to S2 stimulus. See below for advice.
#' @param s Ratio of standard deviations for the S1 and S2 stimulus. Defaults to
#'   1.
#' @param add_constant Adds a small constant to the data (1/number of possible
#'   responses) to account for 0 or 1 values. Defaults to TRUE for ease of use
#'   across multiple datasets.
#'
#' @author Maniscalco & Lau. Ported to R by Matt Craddock,
#'   \email{matt@mattcraddock.com}
#' @references Maniscalco, B., & Lau, H. (2012). A signal detection theoretic
#'   approach for estimating metacognitive sensitivity from confidence ratings.
#'   Consciousness and Cognition. http://dx.doi.org/10.1016/j.concog.2011.09.021
#' @import dplyr
#' @import tidyr
#' @importFrom stats optim pnorm qnorm
#' @examples
#' nR_S1 <- c(60,30,10,7,4,1)
#' nR_S2 <- c(4,6,11,13,23,61)
#' fit_meta_d_MLE(nR_S1, nR_S2)
#' @export
#'

fit_meta_d_MLE <- function(nR_S1, nR_S2, s = 1, add_constant = TRUE) {

  if (add_constant) {
    nR_S1 <- nR_S1 + (1 / length(nR_S1))
    nR_S2 <- nR_S2 + (1 / length(nR_S2))
  }

  n_ratings <- length(nR_S1) / 2
  n_criteria <- 2 * n_ratings - 1

  A <- matrix(0, nrow = n_criteria - 2, ncol = n_criteria)

  for (i in 2:n_criteria - 1) {
    A[i - 1, i:(i + 1)] <- c(1, -1)
  }

  b <- rep(-1e-5, n_criteria - 2)

  ## set up initial guess at parameter values
  ratingHR  <- NULL
  ratingFAR <- NULL

  for (c in 2:(n_ratings * 2)) {
    ratingHR[c - 1] <- sum(nR_S2[c:length(nR_S2)]) / sum(nR_S2)
    ratingFAR[c - 1] <- sum(nR_S1[c:length(nR_S1)]) / sum(nR_S1)
  }

  t1_index <- n_ratings
  t2_index <- setdiff(1:(2 * n_ratings - 1), t1_index)

  d1 <- (1 / s) * qnorm(ratingHR[t1_index]) - qnorm(ratingFAR[t1_index])
  meta_d1 <- d1

  c1 <- (-1 / (1 + s)) * (qnorm(ratingHR) + qnorm(ratingFAR) )
  t1c1 <- c1[t1_index]
  t2c1 <- c1[t2_index]

  guess <- c(meta_d1, t2c1 - (meta_d1 * t1c1 / d1))

  params <- list("n_ratings" = n_ratings,
                 "d1" = d1,
                 "nR_S1" = nR_S1,
                 "nR_S2" = nR_S2,
                 "t1c1" = t1c1,
                 "s" = s)
  #Note that I suppress warnings from the optimizer about NaNs generated while
  #calculating the log-likelihood. They don't effect the final estimates anyway,
  #so the warnings are just annoying

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

  C_nR_rS2 <- nR_S2[(n_ratings + 1):length(nR_S2)]
  C_nR_rS1 <- nR_S1[seq(n_ratings, 1, -1)]
  obs_FAR2_rS2 <- matrix(0, nrow = n_ratings - 1)
  obs_HR2_rS2 <- matrix(0, nrow = n_ratings - 1)
  obs_FAR2_rS1 <- matrix(0, nrow = n_ratings - 1)
  obs_HR2_rS1 <- matrix(0, nrow = n_ratings - 1)

  for (i in 2:n_ratings) {
    obs_FAR2_rS2[i - 1] <- sum(I_nR_rS2[i:length(I_nR_rS2)] ) / sum(I_nR_rS2)
    obs_HR2_rS2[i - 1]  <- sum(C_nR_rS2[i:length(I_nR_rS2)] ) / sum(C_nR_rS2)

    obs_FAR2_rS1[i - 1] <- sum(I_nR_rS1[i:length(I_nR_rS2)] ) / sum(I_nR_rS1)
    obs_HR2_rS1[i - 1]  <- sum(C_nR_rS1[i:length(I_nR_rS2)] ) / sum(C_nR_rS1)
  }

  S1mu <- -meta_d1 / 2
  S1sd <- 1
  S2mu <-  meta_d1 / 2
  S2sd <- S1sd / s
  mt1c1 <- (meta_d1 * t1c1 / d1);

  C_area_rS2 <- 1 - pnorm(mt1c1, S2mu, S2sd)
  I_area_rS2 <- 1 - pnorm(mt1c1, S1mu, S1sd)
  C_area_rS1 <- pnorm(mt1c1, S1mu, S1sd)
  I_area_rS1 <- pnorm(mt1c1, S2mu, S2sd)

  est_FAR2_rS2 <- NULL
  est_FAR2_rS1 <- NULL
  est_HR2_rS2 <- NULL
  est_HR2_rS1 <- NULL

  for (i in 1:(n_ratings - 1)) {

    t2c1_lower <- t2c1[n_ratings - i]
    t2c1_upper <- t2c1[n_ratings - 1 + i]

    I_FAR_area_rS2 <- 1 - pnorm(t2c1_upper, S1mu, S1sd)
    C_HR_area_rS2  <- 1 - pnorm(t2c1_upper, S2mu, S2sd)

    I_FAR_area_rS1 <- pnorm(t2c1_lower, S2mu, S2sd)
    C_HR_area_rS1  <- pnorm(t2c1_lower, S1mu, S1sd)

    est_FAR2_rS2[i] <- I_FAR_area_rS2 / I_area_rS2
    est_HR2_rS2[i]  <- C_HR_area_rS2 / C_area_rS2

    est_FAR2_rS1[i] <- I_FAR_area_rS1 / I_area_rS1
    est_HR2_rS1[i]  <- C_HR_area_rS1 / C_area_rS1
  }

  da <- sqrt(2 / (1 + s ^ 2)) * s * d1
  meta_da <- sqrt(2 / (1 + s ^ 2)) * s * meta_d1
  mt1c1 <- (meta_d1 * t1c1 / d1)
  t2ca <- (sqrt(2) * s / sqrt(1 + s ^ 2)) * t2c1

  #Note that the S1units are not currently in the output - may change it to
  #include them

  S1units <- data.frame(d1 = d1,
                        meta_d1 = meta_d1,
                        s = s,
                        meta_c1 = mt1c1,
                        t2c1_rS1 = t2c1[1:n_ratings - 1],
                        t2c1_rS2 = t2c1[n_ratings:length(t2c1)])

  fit <- data.frame(da = da,
                    s = s,
                    meta_da = meta_da,
                    M_diff = meta_da - da,
                    M_ratio = meta_da / da,
                    meta_ca = (sqrt(2) * s / sqrt(1 + s ^ 2) * mt1c1),
                    t2ca_rS1 = t2ca[1:n_ratings - 1],
                    t2ca_rS2 = t2ca[n_ratings:length(t2ca)],
                    logL = logL,
                    est_HR2_rS1 = est_HR2_rS1,
                    est_HR2_rS2 = est_HR2_rS2,
                    est_FAR2_rS1 = est_FAR2_rS1,
                    est_FAR2_rS2 = est_FAR2_rS2,
                    obs_HR2_rS1 = obs_HR2_rS1,
                    obs_HR2_rS2 = obs_HR2_rS2,
                    obs_FAR2_rS1 = obs_FAR2_rS1,
                    obs_FAR2_rS2 = obs_FAR2_rS2,
                    t1c1 = t1c1
  )
  return(fit)
}


#' Likelihood function for fitting meta-d.
#'
#' This is a likelihood function for use in MLE estimation, and shouldn't be
#' called directly.
#'
#' @param x Starting guess for parameter values
#' @param parameters Various parameters such as the number of ratings, type 1
#'   d-prime etc.
#' @author Maniscalco and Lau. Ported to R by Matt Craddock,
#'   \email{matt@@mattcraddock.com}
#'

fit_meta_d_logL <- function(x, parameters) {
  meta_d1 <- x[1]
  t2c1 <- x[2:length(x)]

  S1mu <- -meta_d1 / 2
  S1sd <- 1
  S2mu <- meta_d1 / 2
  S2sd <- S1sd / parameters$s

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
  C_area_rS1 <- pnorm(t1c1, S1mu, S1sd)
  I_area_rS1 <- pnorm(t1c1, S2mu, S2sd)

  C_area_rS2 <- 1 - pnorm(t1c1, S2mu, S2sd)
  I_area_rS2 <- 1 - pnorm(t1c1, S1mu, S1sd)

  t2c1x <- c(-Inf, t2c1[1:parameters$n_ratings - 1], t1c1,
             t2c1[parameters$n_ratings:length(t2c1)], Inf)
  prC_rS1 <- matrix(0, ncol = parameters$n_ratings)
  prI_rS1 <- matrix(0, ncol = parameters$n_ratings)
  prC_rS2 <- matrix(0, ncol = parameters$n_ratings)
  prI_rS2 <- matrix(0, ncol = parameters$n_ratings)

  for (i in 1:parameters$n_ratings) {
    prC_rS1[i] <- (pnorm(t2c1x[i + 1], S1mu, S1sd) -
                     pnorm(t2c1x[i], S1mu, S1sd)) / C_area_rS1
    prI_rS1[i] <- (pnorm(t2c1x[i + 1], S2mu, S2sd) -
                     pnorm(t2c1x[i], S2mu, S2sd)) / I_area_rS1

    prC_rS2[i] <-
      ((1 - pnorm(t2c1x[parameters$n_ratings + i], S2mu, S2sd)) -
         (1 - pnorm(t2c1x[parameters$n_ratings + i + 1], S2mu, S2sd))) /
      C_area_rS2
    prI_rS2[i] <-
      ((1 - pnorm(t2c1x[parameters$n_ratings + i], S1mu, S1sd)) -
         (1 - pnorm(t2c1x[parameters$n_ratings + i + 1], S1mu, S1sd))) /
      I_area_rS2
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


#'Function for calculating meta-d' by minimizing SSE
#'
#'Provides a type-2 SDT analysis of data from a typical experiment in which
#'observers discriminate between two response alternatives and provide ratings
#'of confidence in their judgements.
#'
#'Where \code{fit_meta_d_MLE} uses Maximum Likelihood Estimation,
#'\code{fit_meta_d_SSE} works by finding the minimum sum of squared errors. As
#'with the MLE method, input is expected as counts for each of two stimulus
#'types.
#'
#'The expected input is two vectors, one for responses to each stimulus,
#'encoding the observers response and confidence. For example, for two stimului
#'labelled A and B, with three confidence ratings, participants could respond to
#'stimulus A as follows: Response: A, rating: 3, count: 60 Response: A, rating:
#'2, count: 30 Response: A, rating: 1, count: 10 Response: B, rating: 1, count:
#'7 Response: B, rating: 2, count: 4 Response: B, rating: 3, count: 1
#'
#'The appropriate vector would be nR_S1 <- c(60,30,10,7,4,1)
#'
#'For stimulus B, we would have the respective vector for responses to stimulus
#'B, eg: Response: A, rating: 3, count: 4 Response: A, rating: 2, count: 6
#'Response: A, rating: 1, count: 11 Response: B, rating: 1, count: 13 Response:
#'B, rating: 2, count: 23 Response: B, rating: 3, count: 61
#'
#'nR_S2 <- c(4,6,11,13,23,61)
#'
#'The output is a dataframe with various metacognitive measures, including
#'m-ratio and meta-d, estimated through minimization of SSE.
#'
#'Currently, multiple rows will be returned when there are more than 2
#'confidence ratings, with some values varying, as they represent cutpoints
#'between confidence ratings, and others simply duplicated.
#'
#'For more details, see Maniscalco & Lau's webpage
#'http://www.columbia.edu/~bsm2105/type2sdt/ Please cite that page and their
#'articles if using this command.
#'
#'@param nR_S1 Responses to S1 stimulus. See below for advice.
#'@param nR_S2 Responses to S2 stimulus. See below for advice.
#'@param s Ratio of standard deviations for the S1 and S2 stimulus. Defaults to
#'  1.
#'@param d_min Minimum bound for d'
#'@param d_max Maximum bound for d'
#'@param d_grain Resolution of grid of possible parameters between the bounds.
#'@param add_constant Adds a small constant to the data (1/number of possible
#'  responses) to account for 0 or 1 values. Defaults to TRUE for ease of use
#'  across multiple datasets.
#'@author Maniscalco & Lau. Ported to R by Matt Craddock
#'  \email{matt@@mattcraddock.com}
#'@import dplyr
#'@importFrom purrr map_dbl map2 map
#' @family [meta_d]
#' @examples
#' nR_S1 <- c(60,30,10,7,4,1)
#' nR_S2 <- c(4,6,11,13,23,61)
#' fit_meta_d_SSE(nR_S1, nR_S2)
#'@export

fit_meta_d_SSE <- function(nR_S1, nR_S2, s = 1, d_min = -5, d_max = 5,
                           d_grain = .01, add_constant = TRUE) {

  if (add_constant) {
    nR_S1 <- nR_S1 + (1 / length(nR_S1))
    nR_S2 <- nR_S2 + (1 / length(nR_S2))
  }

  n_ratings <- length(nR_S1) / 2
  n_criteria <- 2 * n_ratings - 1

  S1_HR <- sum(nR_S2[(n_ratings + 1):(n_ratings * 2)]) / sum(nR_S2)
  S1_FA <- sum(nR_S1[(n_ratings + 1):(n_ratings * 2)]) / sum(nR_S1)

  d_1 <- (1 / s) * qnorm(S1_HR) - qnorm(S1_FA)
  c_1 <- (-1 / (1 + s)) * (qnorm(S1_HR) + qnorm(S1_FA))
  c_prime <- c_1 / d_1

  obs_HR2_rS1 <- NULL
  obs_FAR2_rS1 <- NULL

  for (i in 1:(n_ratings - 1)) {
    obs_HR2_rS1[i] <- sum(nR_S1[1:i]) / sum(nR_S1[1:n_ratings])
    obs_FAR2_rS1[i] <- sum(nR_S2[1:i]) / sum(nR_S2[1:n_ratings])
  }

  obs_HR2_rS2 <- NULL
  obs_FAR2_rS2 <- NULL
  for (i in (n_ratings + 2):(2 * n_ratings)) {
    obs_HR2_rS2[(i - (n_ratings + 2) + 1)] <-
      sum(nR_S2[i:(n_ratings * 2)]) /
      sum(nR_S2[(n_ratings + 1):(n_ratings * 2)])

    obs_FAR2_rS2[(i - (n_ratings + 2) + 1)] <-
      sum(nR_S1[i:(n_ratings * 2)]) /
      sum(nR_S1[(n_ratings + 1):(n_ratings * 2)])
  }

  d_grid <- seq(d_min,
                d_max,
                by = d_grain)
  c_grid <- c_prime * d_grid

  S1mu <- -d_grid / 2
  S2mu <- d_grid / 2
  S1sd <- 1
  S2sd <- 1 / s

  bounds <- 5 * max(S1sd,
                    S2sd)
  SSEmin <- Inf

  param_space <- purrr::map2(S1mu,
                            S2mu,
                      ~seq(.x - bounds,
                           .y + bounds,
                           by = .001)) # param space

  min_z <- purrr::map2(param_space,
                       c_grid,
                       ~min(abs(.x - .y)))
  c_ind <- purrr::map2(param_space,
                       c_grid,
                      ~which.min(abs(.x - .y)))

  HRs <- purrr::map2(param_space, S2mu,
                     ~ 1 - (pnorm(.x, .y, S2sd)))
  FARs <- purrr::map2(param_space, S1mu,
                      ~ 1 - (pnorm(.x, .y, S1sd)))

  # fit type 2 data for S1 responses
  est_HR2s_rS1 <- purrr::map2(FARs,
                              c_ind,
                              ~(1 - .x[1:.y]) / (1 - .x[.y]))
  est_FAR2s_rS1 <- purrr::map2(HRs,
                               c_ind,
                               ~(1 - .x[1:.y]) / (1 - .x[.y]))

  SSE <- NULL
  SSE_rS1 <- matrix(data = NaN,
                    nrow = length(d_grid),
                    ncol = n_ratings - 1)
  rS1_ind <- matrix(data = NaN,
                    nrow = length(d_grid),
                    ncol = n_ratings - 1)

  for (n in (1:(n_ratings - 1))) {
    SSE <-
      purrr::map2(est_HR2s_rS1,
                  est_FAR2s_rS1,
                  ~ (.x - obs_HR2_rS1[[n]]) ^ 2 + (.y - obs_FAR2_rS1[[n]]) ^ 2)
    SSE_rS1[, n] <- purrr::map_dbl(SSE,
                                   min)
    inds <- unlist(purrr::map(SSE,
                              which.min))
    rS1_ind[1:length(inds), n] <- inds
  }

  # fit type 2 data for S2 responses
  est_HR2s_rS2 <- purrr::map2(HRs,
                              c_ind,
                              ~ .x[.y:length(.x)] / .x[.y])
  est_FAR2s_rS2 <- purrr::map2(FARs,
                               c_ind,
                               ~ .x[.y:length(.x)] / .x[.y])

  SSE_rS2 <- matrix(data = NaN,
                    nrow = length(d_grid),
                    ncol = n_ratings - 1)
  rS2_ind <- matrix(data = NaN,
                    nrow = length(d_grid),
                    ncol = n_ratings - 1)
  for (n in (1:(n_ratings - 1))) {
    SSE <- purrr::map2(est_HR2s_rS2,
                       est_FAR2s_rS2,
                       ~(.x - obs_HR2_rS2[[n]]) ^ 2 + (.y - obs_FAR2_rS2[[n]]) ^ 2)
    SSE_rS2[, n] <- purrr::map_dbl(SSE, min)
    inds <- unlist(purrr::map(SSE,
                              which.min))
    rS2_ind[1:length(inds), n] <- inds
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
                    obs_FAR2_rS2 = obs_FAR2_rS2,
                    t1c1 = c_1)
  return(out)
}


#' Meta-d' balance calculation
#'
#' The expected input is two vectors, one for responses to each stimulus,
#' encoding the observers response and confidence. For example, for two stimului
#' labelled A and B, with three confidence ratings, participants could respond
#' to stimulus A as follows: Response: A, rating: 3, count: 60 Response: A,
#' rating: 2, count: 30 Response: A, rating: 1, count: 10 Response: B, rating:
#' 1, count: 7 Response: B, rating: 2, count: 4 Response: B, rating: 3, count: 1
#'
#' The appropriate vector would be nR_S1 <- c(60,30,10,7,4,1)
#'
#' For stimulus B, we would have the respective vector for responses to stimulus
#' B, eg: Response: A, rating: 3, count: 4 Response: A, rating: 2, count: 6
#' Response: A, rating: 1, count: 11 Response: B, rating: 1, count: 13 Response:
#' B, rating: 2, count: 23 Response: B, rating: 3, count: 61
#'
#' nR_S2 <- c(4,6,11,13,23,61)

#' @param nR_S1 Responses to S1 stimulus. See below for advice.
#' @param nR_S2 Responses to S2 stimulus. See below for advice.
#' @param s Ratio of standard deviations for the S1 and S2 stimulus. Defaults to
#'   1.
#' @param add_constant Add a small constant to all cells to adjust for boundary
#'   issues, and for consistency with the use of this method with other meta-d
#'   measures. Note: default for this is FALSE.
#'
#' @author Adam Barrett. Ported to R by Matt Craddock
#'   \email{matt@mattcraddock.com}
#' @references Barrett, Dienes, & Seth (2013). Measures of metacognition on
#'   signal-detection theoretic models. Psychol Methods, 18.
#'   http://dx.doi.org/10.1037/a0033268
#' @importFrom nleqslv nleqslv
#' @family [meta_d]
#' @examples
#' nR_S1 <- c(60,30,10,7,4,1)
#' nR_S2 <- c(4,6,11,13,23,61)
#' fit_meta_d_bal(nR_S1, nR_S2)
#' @export
#'

fit_meta_d_bal <- function(nR_S1,
                           nR_S2,
                           s = 1,
                           add_constant = FALSE) {

  if (add_constant) {
    nR_S1 <- nR_S1 + (1 / length(nR_S1))
    nR_S2 <- nR_S2 + (1 / length(nR_S2))
  }

  n_ratings <- length(nR_S1) / 2

  S1_HR <- sum(nR_S1[1:n_ratings]) / sum(nR_S1)
  S1_FA <- sum(nR_S2[1:n_ratings]) / sum(nR_S2)

  d_prime <- (1 / s) * qnorm(S1_HR) - qnorm(S1_FA)
  t1c1 <- (-1 / (1 + s)) * (qnorm(S1_HR) + qnorm(S1_FA))

  Hp <- nR_S1[[1]] / sum(nR_S1[1:n_ratings])
  Hm <- nR_S2[[n_ratings * 2]] / sum(nR_S2[(n_ratings + 1):(n_ratings * 2)])
  Fm <- nR_S1[[n_ratings * 2]] / sum(nR_S1[(n_ratings + 1):(n_ratings * 2)])
  Fp <- nR_S2[[1]] / sum(nR_S2[1:n_ratings])

  theta <- -qnorm(S1_FA)
  theta_prime <- theta / d_prime
  x0 <- c(theta, d_prime)

  ep <- nleqslv::nleqslv(x0, fit_metad_plus, method = "Newton", th = theta_prime,
                hp = Hp, fp = Fp, global = "pwldog",
                control = list(delta = "cauchy", ftol = 1e-06), xscalm = "auto")

  meta_d_plus <- ep$x[[2]]
  x0 <- c(theta_prime, d_prime)

  em <- nleqslv(x0, fit_metad_minus, method = "Newton", th = theta_prime,
                hm = Hm, fm = Fm, global = "pwldog",
                control = list(delta = "cauchy", ftol = 1e-06), xscalm = "auto")
  meta_d_neg <- em$x[[2]]

  htp <- 1 - pnorm(theta_prime * meta_d_plus, meta_d_plus, s)
  ftp <- 1 - pnorm(theta_prime * meta_d_plus, 0, 1)
  htm <- 1 - pnorm(theta_prime * meta_d_neg, meta_d_neg, s)
  ftm <- 1 - pnorm(theta_prime * meta_d_neg, 0, 1)

  if (any(c(htp, ftp, htm, ftm) > .95) | any(c(htp, ftp, htm, ftm) < .05)) {
    stable <- 0
  } else {
    stable <- 1
  }

  r <- (S1_HR + S1_FA) / 2
  meta_d_bal <- r * meta_d_plus + (1 - r) * meta_d_neg

  data.frame(d_prime = d_prime,
             meta_d_bal = meta_d_bal,
             meta_d_ratio = meta_d_bal / d_prime,
             meta_d_diff = meta_d_bal - d_prime,
             stable = stable,
             htp = htp,
             ftp = ftp,
             htm = htm,
             ftm = ftm,
             t1c1 = t1c1)
}

#' Internal function for fitting meta_d_plus
#'
#' @param x0 Starting parameters.
#' @param th Theta.
#' @param hp Hit rate for positive responses
#' @param fp False positive rate for positive responses.
#' @keywords internal

fit_metad_plus <- function(x0, th, hp, fp) {
  y1 <- (1 - pnorm(x0[[1]], x0[[2]], 1)) /
    (1 - pnorm(th * x0[[2]], x0[[2]], 1)) - hp
  y2 <- (1 - pnorm(x0[[1]], 0, 1)) / (1 - pnorm(th * x0[[2]], 0, 1)) - fp
  c(y1, y2)
}


#' Internal function for fitting meta_d_minus
#'
#' @param x0 Starting parameters.
#' @param th Theta
#' @param hm Hit rate for negative responses
#' @param fm False positive rate for negative responses
#' @keywords internal
fit_metad_minus <- function(x0, th, hm, fm) {
  y1 <- pnorm(x0[[1]], 0, 1) / pnorm(th * x0[[2]], 0, 1) - hm
  y2 <- pnorm(x0[[1]], x0[[2]], 1) / pnorm(th * x0[[2]], x0[[2]], 1) - fm
  c(y1, y2)
}
