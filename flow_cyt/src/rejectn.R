library('MASS')

get_sample <- function(mu = rep(0,4), Sig = diag(4)) {
  smpls <- mvrnorm(1, mu, Sig)
  acc   <- isconvex(smpls)
  while(!acc) {
    smp <- mvrnorm(1, mu, Sig)
    if(isconvex(smp)) {
      smpls <- rbind(smpls, smp)
    } else smpls <- rbind(smpls, smp)
    acc   <- isconvex(smp)
  }
  return(smpls)
}


get_sample_smplx <- function(mu, Sig) {
  acc <- F
  while(acc) {
    smp <- mvrnorm(blk, mu, Sig)
    tst <- condn(smp)
    if(any(tst)) {
      smpls <- rbind(smpls, smp[1:which(tst),])
    } else smpls <- rbind(smpls, smp)
  }
}


get_sample_smplx <- function(mu, Sig) {
  acc <- F
  while(acc) {
    smp <- mvrnorm(blk, mu, Sig)
    tst <- condn(smp)
    if(any(tst)) {
      smpls <- rbind(smpls, smp[1:which(tst),])
    } else smpls <- rbind(smpls, smp)
  }
}

condn <- function(inp) {
  (apply(inp,1,min) > 0) & (rowSums(jnk) < 1)
}

isconvex <- function(inp) {
  if(length(inp)==2) return(T)

  ln   <- length(inp)
  dif  <- inp[2:ln] - inp[1:(ln-1)]
  dif2 <- dif[2:(ln-1)] - dif[1:(ln-2)]
  all(dif2>0)
}
