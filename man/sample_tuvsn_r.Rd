% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samp_cases.R
\name{sample_tuvsn_r}
\alias{sample_tuvsn_r}
\title{Full Rejection sampling steps}
\usage{
sample_tuvsn_r(a, b)
}
\arguments{
\item{a}{Lower bound of truncation}

\item{b}{Upper bound of truncattion}
}
\value{
\item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
}
\description{
Generate a sample from a truncated univariate standard normal distribution,
    TN(0, 1; a, b).
}
\references{
Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
    truncated multivariate normal and student-t distributions subject to
    linear inequality constraints. Journal of Statistical Theory and
    Practice, 9(4), 712-732.
}

