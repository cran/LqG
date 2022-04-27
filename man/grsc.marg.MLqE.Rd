\name{grsc.marg.MLqE}
\alias{grsc.marg.MLqE}
\title{Group Screening based on marginal Maximum Lq-likelihood Estimation}
\description{
 Group screening by ranking utility of each group. The group effect is defined based on the cumulation of the maximum Lq-likelihood estimate of the regression using only one predictor each time within the group.
}
\usage{
grsc.marg.MLqE(
X,
Y,
n = dim(X)[1],
p = dim(X)[2],
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
\item{p}{A value denoting the dimension of predictors}
\item{q}{A value of distortion parameter of Lq function, default to \code{0.9}.}
\item{m}{A number of the predictor groups}
\item{group}{A vector of consecutive integers describing the grouping of the coefficients (see example below).}
\item{eps}{The iteration coverage criterion, default to \code{1e-06}.}
\item{d}{A value of the number of groups retained after screening, default to \code{n/log(n)}}.}
\details{ grsc.marg.MLqE obtains the group effect of each group for subsequential group screening, based on the cumulative marginal MLqE coefficients within the group. It can work when both the correlation within groups and between groups are small. If group size equals to 1, individual screening is conducted. }
\value{
The \code{grsc.marg.MLqE} returns a list containing the following components:
\item{beta.group}{The vector of utility of each group, which is the criterion for the variable screening procedure.}
\item{group.screened}{The vector of integers denoting the screened groups.}}
\examples{
# This is an example of grsc.marg.MLqE with simulated data
data(LqG_SimuData)
X = LqG_SimuData$X
Y = LqG_SimuData$Y
n = dim(X)[1]
p = dim(X)[2]
m = 200
groups = rep(1:(p/5), each = 5)
result <- grsc.marg.MLqE(X = X,
                         Y = Y,
                         n = n,
                         p = p,
                         q = 0.9,
                         m = m,
                         group = groups,
                         eps = 1e-06,
                         d = 15)
result$beta.group
result$group.screened
}


