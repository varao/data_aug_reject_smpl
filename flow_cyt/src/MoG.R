library('MCMCpack')

mog <- function(ip_data, K) {

  num_iter <- 1000
  alp      <- .1

  dm  <- dim(ip_data)[1]
  N   <- dim(ip_data)[2]

  niw_p <- list(mu0 = rep(0, dm), lam = .1, Phi = diag(rep(dm)),  nu = dm+1)

  wt  <- rep(0, K)

  mu  <- list()
  cv  <- list()
# mu[[1]] <- c(2,2)
# mu[[2]] <- c(-3,-3)

# cv[[1]] <- diag(2)
# cv[[2]] <- diag(2)

  ##################
  # Initialize
  z  <- sample.int(K, N, replace=T)

  for(iter in 1:num_iter) {

    print(iter)
    params <- updt_params(K, alp, ip_data, dm, niw_p, z)

    z <- updt_assgn(K, ip_data, N, dm, params) 
  }
  print(mu)
  print(cv)
  return(list(z=z,mu=mu,cv=cv, wt = wt))
}

updt_params <- function(K, alp, ip_data, dm, niw_p, z) { 

  mu  <- list()
  cv  <- list()

  count  <- rep(0, K)
  for(i in 1:K) count[i] <- sum(z==i)
  wt <- rgamma(K,alp+count, 1)
  wt <- wt / sum(wt)

# count  <- tabulate(z,K)
# cum    <- rev(cumsum(rev(count)))
# bt     <- rbeta(K-1, 1+count, alp+c(cum[2:K]))
# bt[K]  <- 1
# cmp_bt <- 1-bt
# wt     <- bt
# for(i in 2:K) 
#   wt[i] <- prod(cmp_bt[1:(i-1)])*bt[i]

  for(i in 1:K) {
    nm   <- count[i]
    X    <- ip_data[,z==i, drop=F]
    if(nm > 0) {
      mn     <- rowMeans(X)
      S      <- (X-mn) %*% t(X-mn)
    } else {
      mn     <- rep(0,dm)
      S      <- diag(rep(0, dm))
    }
    SS   <- list(n=nm, mn=mn, S=S)
    rslt <- niw.post(niw_p$mu0, niw_p$lam, niw_p$Phi, niw_p$nu, SS) 
    mu[[i]] <- rslt$mu
    cv[[i]] <- rslt$cv
  }
  list(wt = wt, mu = mu, cv = cv)
}

updt_assgn <- function(K, ip_data, N, dm, params) {

  wt <- params$wt
  mu <- params$mu
  cv <- params$cv

  lik <- matrix(rep(0,N*K), nrow=K)
  for(i in 1:K) {
    #lik[,i] <- dnorm(ip_data, mu[[i]], cv[[i]], log=T)
    #print(dim(((ip_off %*% solve(cv[[i]], t(ip_off))))))
    ip_off  <- ip_data - mu[[i]]
    rot <- solve(cv[[i]], ip_off)
    lik[i,] <- -0.5*colSums(ip_off * rot) - 0.5 * determinant(cv[[i]],log=T)$modulus
  }
  lik <- lik + log(wt)
  lik <- lik - apply(lik,2,max)
  lik <- t(t(lik) - log(colSums(exp(lik))))
  lik <- exp(lik)
  lik <- apply(lik,2,cumsum)
  z <- rowSums(runif(N) > t(lik)) + 1L
}


gen_mog <- function(K = 2, params,  N= 100) {
  wt  <- params$wt
  mu  <- params$mu
  cv  <- params$cv

  dm <- length(mu[[1]])

  X <- matrix(rep(0,N*dm), ncol=dm)  # Will return transpose
  z <- sample(K, N, replace=T, prob=wt)

  for(i in 1:K) {
    indx  <- z == i
    if(sum(indx) > 0) X[indx,] <- rmvnorm(sum(indx), mu[[i]], cv[[i]]) 
  }
  list(x=t(X),z=z)
}
