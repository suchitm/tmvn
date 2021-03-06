% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling_methods.R
\name{exp_rej_r}
\alias{exp_rej_r}
\title{TUVN exponential rejection sampling}
\usage{
exp_rej_r(a, b = Inf, lambda = NULL)
}
\arguments{
\item{a}{Lower bound of truncation}

\item{b}{Upper bound of truncation; default is infinity}

\item{lambda}{rate parameter for the exponential distribution. Optimal value is
calculated by default. Refer to Lemma 2.1.1.2 in Li and Ghosh (2015).}
}
\value{
\item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
}
\description{
TUVN exponential rejection sampling
}
\references{
Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
    truncated multivariate normal and student-t distributions subject to
    linear inequality constraints. Journal of Statistical Theory and
    Practice, 9(4), 712-732.
}

