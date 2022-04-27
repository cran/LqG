\name{MLqE.est}
\alias{MLqE.est}
\title{Maximum Lq-likelihood Estimation}
\description{
 The iterative algorithm for MLqE of coefficients of regression using each group of variables.
}
\usage{
MLqE.est(
X,
Y,
q = 0.9,
eps = 1e-06
)
}
\arguments{
\item{X}{The matrix of the predictor group.}
\item{Y}{The vector of response.}
\item{q}{The value of distortion parameter of Lq function, default to \code{0.9}.}
\item{eps}{The iteration coverage criterion, default to \code{1e-06}.}
}
\details{ The estimating equation of MLqE is a weighted version of that of the classical maximum likelihood estimation (MLE) where
the distortion parameter q determines the similarity between the Lq function and the log function. When q = 1, MLqE is equivalent to MLE. The closer q is to 1, the more sensitive the MLqE is to outliers. As for the selection of q, there is presently no general method. However, MLqE is generally less sensitive to data contamination than MLE (to different degrees) when q is smaller than 1. Here, the default value of q is 0.9. Distortion parameter q can also be determined according to sample size n, choices of \eqn{q_n} with \eqn{|1-q_n|} between \eqn{\frac{1}{n}} and \eqn{\frac{1}{\sqrt{n}}} usually improves over the MLE.
}
\value{
The \code{MLqE.est} returns a list containing the following components:
\item{t}{The integer specifying the number of the total iterations in the algorithm.}
\item{beta_hat}{The vector of estimated coefficients.}
\item{sigma_hat}{The value of the estimated variance.}
\item{OMEGA_hat}{The matrix of the estimated weight.}}
\examples{
# This is an example of grsc.marg.MLqE with simulated data
data(LqG_SimuData)
X = LqG_SimuData$X
Y = LqG_SimuData$Y
n = dim(X)[1]
p = dim(X)[2]
m = 200
groups = rep(1:( dim(X)[2] / 5), each = 5)
Xb = X[ , which( groups == 1)]
result = MLqE.est(Xb,
                  Y,
                  q = 0.9,
                  eps = 1e-06)
result$beta_hat
result$sigma_hat
result$OMEGA_hat
result$t
}
