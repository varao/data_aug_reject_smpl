inf_gvhd <- function(K) {
  axs <- 1
  sc <- 1024/axs

# ip_data <- t(gvhd10[,c(1,2,3)])
# ip_data <- ip_data[,!(ip_data[2,]==1024)]
# ip_data <- ip_data[,!(ip_data[3,]==1024)]/sc
# nr      <- sqrt(ip_data[1,]*ip_data[1,] + ip_data[2,]*ip_data[2,]) < 0.75

# print(sum(nr))
# nr      <- nr & (runif(length(nr)) > 0.00)
# print(sum(nr))
# ip_data <- ip_data[,!nr]
# ip_data <- t(kmeans(t(ip_data),1000)$centers)/sc  #, 2000, nstart=1, iter.max=10)$centers)

  ip_data <- t(GvHD.pos)/sc
  dm      <- dim(ip_data)[1]
  N       <- dim(ip_data)[2]
  print(N)
  print(dm)

  l1 <- 0/sc
  u  <- 1023/sc

  genconst <- function(u1,u2,l1,l2) { 
    in_set <- function(X) {
      X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2 & X[3,] > l1 & X[3,] < u1 & X[4,] > l2 & X[4,] < u2 
    }
  }
  constr  <- genconst(u,u,l1,0)
  constr2 <- genconst(1.1*u,1.1*u,-u/10,-u/10)


  num_iter <- 1000
  alp      <- 1
  niw_p <- list(mu0 = rep(axs/2, dm), lam = 1, Phi = diag(rep(1,dm)),  nu = dm+1)
  cnt_rej <- rep(0, num_iter)

  z  <- sample.int(K, N, replace=T)
  z  <- kmeans(t(ip_data), 5, nstart=10, iter.max=50)$cluster
  params  <- updt_params(K, alp, ip_data, dm, niw_p, z)
 # tmp  <- kmeans(t(ip_data), 100, nstart=10, iter.max=50)$centers
 # params$mu  <- kmeans(tmp, K, nstart=10, iter.max=50)$centers
 #params  <- updt_params(K, alp, ip_data, dm, niw_p, z)
 #for(i in 1:K) {
 #  #params$mu[[i]] <- runif(dm)
 #  params$cv[[i]] <- diag(rep(1, dm))
 #}

  for(iter in 1:20) {
    if(iter%%10 == 0) print(iter)
    z       <- updt_assgn(K, ip_data, N, dm, params)
    params  <- updt_params(K, alp, ip_data, dm, niw_p, z)
    #for(i in 1:K) params$cv[[i]] <- diag(rep(.02, dm))
  }

  rslt  <- rep(list(params), num_iter)
  for(iter in 1:num_iter) {
    if(iter%%10 == 0) print(iter)
    z          <- updt_assgn(K, ip_data, N, dm, params)
#    rej_data   <- impute_rej(K, params, N, constr, constr2)
#   params     <- updt_params(K, alp, cbind(ip_data, rej_data$rejs), dm, niw_p, c(z,rej_data$z))
    for(i in 1:K) {
      pm <- updt_params_component(ip_data[,z==i, drop=F], params$mu[[i]], params$cv[[i]], niw_p, constr)
      params$mu[[i]] <- pm$mu
      params$cv[[i]] <- pm$cv
      cnt_rej[iter]  <- cnt_rej[iter] + pm$cnt_rej
    }

    count  <- rep(0, K)
    for(i in 1:K) count[i] <- sum(z==i)
    cum    <- rev(cumsum(rev(count)))
    bt     <- rbeta(K-1, 1+count, alp+c(cum[2:K]))
    bt[K]  <- 1
    cmp_bt <- 1-bt
    wt     <- bt
    for(i in 2:K) 
      wt[i] <- prod(cmp_bt[1:(i-1)])*bt[i]
#   wt <- rgamma(K,alp+count, 1)
#   wt <- wt / sum(wt)
    params$wt <- wt
    #params     <- updt_params(K, alp, ip_data , dm, niw_p, z)
    #for(i in 1:K) params$cv[[i]] <- diag(rep(.02, dm))
    rslt[[iter]] <- params
#    rslt[iter] <- list(data = ip_data, params=params, rej = as.data.frame(rej_data), z = z)
  }
  # Fix the length
  return(list(rslt=rslt, data=ip_data, cnt_rej = cnt_rej))
}

impute_rej <- function(K, params, N, constr, constr2) {
  cum_acc <- 0
  rej_acc <- 0
  rejs   <- matrix(rep(0, 0), nrow = length(params$mu[[1]]))
  z      <- c()
  while(1) {
    smpls  <- gen_mog(K, params, N)
#   vl2     <- constr2(smpls$x)
#   nm     <- sum(vl2)
#   smpls$x <- smpls$x[,vl2]
#   smpls$z <- smpls$z[vl2]
    vl     <- constr(smpls$x)
    num_a  <- sum(vl)
    num_r  <- N - num_a
    rejs   <- cbind(rejs, smpls$x[,!vl])
    z      <- c(z, smpls$z[!vl])
    cum_acc  <- cum_acc + num_a
    rej_acc  <- rej_acc + num_r
    
   
    if(cum_acc >= N) {return(list(rejs=rejs, z=z)) }
  }
}

impute_gauss <- function(mu, cv , N, constr, constr2) {
  cum_acc <- 0
  rej_acc <- 0
  rejs   <- matrix(rep(0, 0), nrow = length(mu))
  z      <- c()
  while(1) {
    smpls <- t(rmvnorm(N, mu, cv))
#   vl2     <- constr2(smpls$x)
#   nm     <- sum(vl2)
#   smpls$x <- smpls$x[,vl2]
#   smpls$z <- smpls$z[vl2]
    vl     <- constr(smpls)
    num_a  <- sum(vl)
    num_r  <- N - num_a
    rejs   <- cbind(rejs, smpls[,!vl])
    cum_acc  <- cum_acc + num_a
    rej_acc  <- rej_acc + num_r
    
    if(cum_acc >= N ) #{print(c(ncol(rejs), cum_acc, N, N*ncol(rejs)/cum_acc)); return(rejs[,floor(N*ncol(rejs)/cum_acc)]) }
        { 
          if(ncol(rejs)==0) return(rejs) else {
            frac_xtr  <- (num_a - (cum_acc - N))/num_a
            len       <- floor(ncol(rejs) - frac_xtr*num_r)
        #   print(c(frac_xtr, len, ncol(rejs)))
           return(rejs[,1:len]) }
        }
        #{return(rejs) }
  }
}

updt_params_component <- function(X, mu, cv, niw_p, constr) { 

  nm   <- ncol(X)
  dm   <- nrow(X)

  if(nm > 0) {
    X <- cbind(X, impute_gauss(mu, cv, nm, constr, constr))
    nm2    <- ncol(X)
    mn     <- rowMeans(X)
    S      <- (X-mn) %*% t(X-mn)
  } else {
    nm2    <- nm
    mn     <- rep(0,dm)
    S      <- diag(rep(0, dm))
  }
  SS   <- list(n=nm2, mn=mn, S=S)
  rslt <- niw.post(niw_p$mu0, niw_p$lam, niw_p$Phi, niw_p$nu, SS) 

  list(mu = rslt$mu, cv = rslt$cv, cnt_rej = nm2-nm)
}
