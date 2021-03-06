% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtmvt.R
\name{rtmvt}
\alias{rtmvt}
\title{Truncated Multivariate Student t Distribution}
\usage{
rtmvt(n, Mean, Sigma, nu, D, lower, upper, init)
}
\arguments{
\item{n}{number of samples to be generated}

\item{Mean}{mean vector}

\item{Sigma}{covariance matrix}

\item{nu}{degress of freedom for the t-distribution}

\item{D}{matrix of linear constraints}

\item{lower}{vector of lower bounds}

\item{upper}{vector of upper bounds}

\item{init}{vector of initial values for the Gibbs sampler. Must satisfy
the linear constraints.}
}
\value{
a matrix of samples with each column being an idependent sample.
}
\description{
This function generates random vectors form the truncates multivariate
    Student t distribution. It uses random sampling from rtmvn since the
    Student t distribution can be represented as a scale-mixture of normals.
}
\examples{
Mean = rep(0,2)
rho = 0.5
Sigma = matrix(c(10,rho,rho,0.1),2,2)
D = matrix(c(1,1,1,-1),2,2)
varp = Sigma[1,1]+Sigma[2,2]+2*Sigma[1,2] # var of the sum
varm = Sigma[1,1]+Sigma[2,2]-2*Sigma[1,2] # var of the diff
sd = c(sqrt(varp),sqrt(varm))
lower = -1.5*sd; upper = 1.5*sd; int = rep(0,2)
nu = 5
n = 20
rtmvt(n,Mean,Sigma,nu,D,lower,upper,int)

}
\references{
Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
    truncated multivariate normal and student-t distributions subject to
    linear inequality constraints. Journal of Statistical Theory and
    Practice, 9(4), 712-732.
}

