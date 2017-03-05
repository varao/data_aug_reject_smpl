library('mixtools')

plot_mog <- function(params, plt) {
  K <- length(params$mu)
  wt <- params$wt / max(params$wt)
  for(i in 1:K) {
#    plt <- plt + geom_point(data=as.data.frame(ellipse(params$mu[[i]][1:2], params$cv[[i]][1:2,1:2])), aes(x=V1,y=V2)) #, alpha=params$wt[[i]])
    plt <- plt + geom_point(data=as.data.frame(ellipse(params$mu[[i]][1:2], params$cv[[i]][1:2,1:2])), aes(x=V1,y=V2), alpha=wt[i])
  }
  return(plt)
}
