## MLqE of coefficients of regression using each group of variables.
MLqE.est <- function(X, Y, q = 0.9, eps = 1e-06){

  beta0 <- solve(t(X)%*%X)%*%t(X)%*%Y
  sigma0 <- crossprod(Y-X%*%beta0)/length(Y)
  t <- 1
  beta_old <- beta0
  sigma_old <- sigma0
  repeat
  {
    omega_hat <- NULL
    for(i in 1:length(Y)){
      omega_hat[i] <- (1/sqrt(2*pi*sigma_old)*exp(-1/(2*sigma_old)*(Y[i]-X[i,]%*%beta_old)^2))^(1-q)
    }
    OMEGA_new <- diag(omega_hat)
    beta_new <- solve(t(X)%*%OMEGA_new%*%X)%*%t(X)%*%OMEGA_new%*%Y
    sigma_new <- sum(omega_hat*(Y-X%*%beta_new)^2)/sum(omega_hat)
    if ((crossprod(beta_new - beta_old) <= eps))
      break
    t <- t + 1
    beta_old <- beta_new
    sigma_old <- sigma_new
  }
  return(list(t=t, beta_hat = beta_new, sigma_hat = sigma_new, OMEGA_hat = OMEGA_new))
}

## Group screening of LqG
grsc.MLqE <- function(X, Y, n = dim(X)[1],
                      q = 0.9, m, group, eps = 1e-06, d = n/log(n)){

  Y <- Y - mean(Y); X <- apply(X, 2, function(c) c-mean(c))
  X <- apply(X, 2, scale)
  beta.group <- NULL
  for(l in 1:m){
    Xb = X[ , which(group==l)]
    betab <- MLqE.est(Xb, Y, q, eps)$beta_hat
    beta_mean <- sqrt(crossprod(betab))/ncol(Xb)
    beta.group <- c(beta.group, beta_mean)
  }
  group.index =  order(beta.group, decreasing = TRUE)[1:floor(d)]

  return(list(beta.group = beta.group,
              group.screened = group.index) )
}

## Group screening of Lq1
grsc.marg.MLqE <- function(X, Y, n = dim(X)[1], p = dim(X)[2],
                           q = 0.9, m, group, eps = 1e-06,  d = n/log(n)){

  Y <- Y - mean(Y); X <- apply(X, 2, function(c) c-mean(c))
  X <- apply(X, 2, scale)
  p = dim(X)[2]
  beta.marginal <- NULL
  for(l in 1:p){
    Xb = as.matrix(X[ ,l])
    betab <- MLqE.est(Xb, Y, q, eps)$beta_hat
    beta.marginal <- c(beta.marginal, betab)
  }
  beta.group <- NULL
  for(l in 1:m){
    beta_mean <- sqrt(crossprod(beta.marginal[which(group==l)]))/length(which(group==l))
    beta.group <- c(beta.group, beta_mean)
  }
  group.index = order(beta.group, decreasing = TRUE)[1:floor(d)]

  return(list(beta.group = beta.group,
              group.screened = group.index) )
}
