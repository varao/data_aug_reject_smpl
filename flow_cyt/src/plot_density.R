plot_density <- function(rslt, num_cls) {
  grd_pt  <- seq(-.5,2,.02)
  grd     <- matrix(rep(0, 2* length(grd_pt)*length(grd_pt)), nrow=2)
  grd[1,] <- rep(grd_pt, length(grd_pt))
  grd[2,] <- rep(grd_pt, each=length(grd_pt))
  grd_wt <- rep(0, ncol(grd))

  dms <- c(1,2)
  for(i in 500:1000) {
    if(i %% 100 == 0) print(i)
    for(k in 1:num_cls) {
      cv <- rslt[[i]]$cv[[k]][dms, dms]
      mn <- rslt[[i]]$mu[[k]][dms]
      
      ip_off  <- grd - mn
      rot <- solve(cv, ip_off)
      tmp <- -0.5*colSums(ip_off * rot) - 0.5 * determinant(cv,log=T)$modulus
      grd_wt <- grd_wt + rslt[[i]]$wt[k] * exp(tmp)
    }
  }
  data.frame(x=grd[1,],y=grd[2,],z=(grd_wt))
}
