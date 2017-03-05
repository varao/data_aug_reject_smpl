library('MCMCpack')

niw <- function(mu0, lam, Phi, nu) {
  cv <- riwish(nu, Phi)
  mu <- mvrnorm(1,mu0, cv / lam)
  return(list(mu=mu,cv=cv))
}

niw.post <- function(mu0, lam, Phi, nu, SS) {
  n     <- SS$n
  mn    <- SS$mn
  S     <- SS$S     # S = t(X - mn) %*% (X - mn)
  mu_n  <- (lam*mu0 + n*mn)/(lam+n)
  Phi_n <- Phi + S + (lam*n)/(lam+n)*outer(mn-mu0,mn-mu0)

  niw(mu_n, lam+n, Phi_n, nu+n)
}

niw.logZ <- function(mu, lam, Phi, nu) {
  logZ <- - 0.5*nu*determinant(Phi,log=T)$modulus + 0.5*nu*dm*log(2) + 0.5*log(2*pi/lam)

  logZ <- logZ + sum(lgamma((dm + 1 - 1:k)/2))
  logZ <- logZ + (dm * (dm - 1)/4)*log(pi)
}

niw.logpx <- function(X, mu0, lam, Phi, nu) {

  mn    <- rowMeans(X)
  n     <- ncol(X)
  S     <- (t(X[,indx]-mn)) %*% (X-mn)
  mu_n  <- (lam*mu0 + n*mn)/(lam+n)
  Phi_n <- Phi + S + (lam*n)/(lam+n)*outer(mn-mu0,mn-mu0)

  niw.logZ(mu_n, lam+n, Phi_n, nu+n)-niw.logZ(mu, lam, Phi, nu, SS)-0.5*nm*log(2*pi)
}

niw.test <- function(dm = 3) {
  nm  <- 100

  mu  <- rep(0, dm)
  lam <- 1
  Phi <- diag(rep(dm))
  nu  <- dm+1

  smp <- niw(mu, lam, Phi, nu)
  print(smp)
  X   <- t(mvrnorm(nm, smp$mu, smp$cv))

  mn <- rowMeans(X)
  SS  <- list(n=nm, mn=mn, S=(X-mn) %*% t(X-mn))
  print(SS$S/SS$n)

  niw.post(mu, lam, Phi, nu, SS) 
}
