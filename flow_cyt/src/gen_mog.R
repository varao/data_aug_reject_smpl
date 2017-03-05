library('MCMCpack')

gen_mog <- function(K = 2, N= 100, dm=2) {
  mu  <- list()
  cv  <- list()

  for(i in 1:K) {
    mu[[i]] <- rnorm(dm)*5
    cv[[i]] <- riwish(dm+1, diag(dm))
  }

X <- matrix(rep(0,N*dm), nrow=dm)
for(i in 1:N) {
  z    <- sample.int(K,1)
  X[,i] <- mvrnorm(1, mu[[z]], cv[[z]]) 
}
return(X)
}
