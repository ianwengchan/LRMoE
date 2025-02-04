% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Predict.R
\name{predict_excess_posterior}
\alias{predict_excess_posterior}
\title{predict_excess_posterior}
\usage{
predict_excess_posterior(
  Y,
  X,
  alpha,
  model,
  limit,
  exposure_past = list(),
  exposure_future = list(),
  exact_Y = F
)
}
\arguments{
\item{Y}{A matrix of responses.}

\item{X}{A matrix of covariates.}

\item{alpha}{A matrix of logit regression coefficients.}

\item{model}{A matrix specifying the expert functions.}

\item{limit}{A matrix specifying the cutoff point.}

\item{exposure_past}{A vector indicating the time exposure (past) of each observation. If nothing is supplied,it is set to 1.0 by default.}

\item{exposure_future}{A vector indicating the time exposure (future) of each observation. If nothing is supplied,it is set to 1.0 by default.}

\item{exact_Y}{Bool variable. indicating if \code{Y} is observed exactly or with censoring and truncation. Default set to be False}
}
\value{
result: A matrix of predicted excess expectation of response, based on posterior probabilities.
}
\description{
Predicts the excess expectation of response, given covariates \code{X},
logit regression coefficients \code{alpha} and a specified \code{model} of expert functions.
}
