crp <- function(X, sz, alp = 1, G0 = c(1,1), num_iter = 1000)
{

  # CRP sampler for DP mixture model with Normal-Inverse-Wishart base-measure
  #
  # Vinayak Rao
  #
  # X is a dataframe with each observation a column 
  # Returns a list of two matrices giving cluster assignment and parameters of each cluster
  #
  # Set the maximum number of clusters in max_cls

  max_cls   <- 100

  num_cls   <- max(X$c)
  num_obs   <- nrow(X)

  ret_c     <- matrix(rep(0,num_obs*num_iter), nrow=num_iter)
  ret_p_bin <- ret_c
  denom     <- num_obs + alp - 1

  n_c       <- rep(0, max_cls)  # Hopefully, 100 > typical num of clusters
  p_vec     <- rep(0, max_cls)
  Z         <- rep(0, max_cls)  # Hopefully, 100 > typical num of clusters

  for(i in 1:num_cls) {
    indx     <- cls[1,] == i
    n_c[i]   <- sum(indx)
    Z[i]     <- niw.logpx(X[,indx], mu0, lam, Phi, nu)
  }

  for(iter in 1:num_iter) {
    if (iter %% 100 == 0) print(iter)
    for(i in 1:num_obs) {

      X_curr      <- [,i]
      c_curr      <- c[iter,i]

      n_c[c_curr] <- n_c[c_curr] - 1
      indx     <- cls[iter,] == c_curr
      Z[c_curr]     <- niw.logpx(X[,indx], mu0, lam, Phi, nu)


      if(n_c[c_curr] == 0) {
        n_c[c_curr]   <- n_c[num_cls]
        n_c[num_cls]  <- 0
        X$c[X$c == num_cls] <- c_curr
        num_cls       <- num_cls - 1    
        # Update param and Z
      }

      p_vec[1:num_cls] <- log(n_c[1:num_cls])
      p_vec[num_cls+1]   <- log(alp)

      c_new        <- sample(num_cls+empt, 1, FALSE, p_vec[indx])

      if(c_new > num_cls) {
        p_bin[num_cls+1] <- p_bin[c_new]
        c_new            <- num_cls+1
        num_cls          <- num_cls + 1
      }
      X$c[i]      <- c_new
      n_c[c_new]  <- n_c[c_new] + 1
   }

  }
  return(list(c=ret_c, p_bin=ret_p_bin))
}
