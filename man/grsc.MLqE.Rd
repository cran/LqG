\name{grsc.MLqE}
\alias{grsc.MLqE}
\title{Group Screening based on Maximum Lq-likelihood Estimation}
\description{
Group screening by ranking utility of each group. The group effect is defined based on the maximum Lq-likelihood estimates of the regression using each group of variables.
}
\usage{
grsc.MLqE(
X,
Y,
n = dim(X)[1],
q = 0.9,
m,
group,
eps = 1e-06,
d = n/log(n)
)
}
\arguments{
\item{X}{A matrix of predictors.}
\item{Y}{A vector of response.}
\item{n}{A value of sample size}
\item{q}{A value of distortion parameter of Lq function, default to \code{0.9}.}
\item{m}{A number of the predictor groups}
\item{group}{A vector of consecutive integers describing the grouping of the coefficients (see example below).}
\item{eps}{The iteration coverage criterion, default to \code{1e-06}.}
\item{d}{A value of the number of groups retained after screening, default to \code{n/log(n)}.} }
\details{ grsc.MLqE obtains the group effect of each group for subsequential group screening, based on the maximum Lq-likelihood estimates of the regression using
each group of variables. By inheriting the advantage of the MLqE in small or moderate sample situations, the method is more robust to heterogeneous data and heavy-tailed distributions. It can work when correlation is mild or large. If group size equals to 1, individual screening is conducted. }
\value{
The \code{grsc.MLqE} returns a list containing the following components:
\item{beta.group}{The vector of utility of each group, which is the criterion for the variable screening procedure.}
\item{group.screened}{The vector of integers denoting the screened groups.}}
\examples{
# This is an example of grsc.MLqE with simulated data
data(LqG_SimuData)
X = LqG_SimuData$X
Y = LqG_SimuData$Y
n = dim(X)[1]
m = 200
groups = rep(1:( dim(X)[2] / 5), each = 5)
result <- grsc.MLqE(X = X,
                    Y = Y,
                    n = n,
                    q = 0.9,
                    m = m,
                    group = groups,
                    eps = 1e-06,
                    d = 15)
result$beta.group
result$group.screened
}

