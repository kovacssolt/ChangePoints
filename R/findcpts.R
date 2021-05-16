#' Changpoint Finder
#'
#' Function that finds changepoints in given data set x by evaluating a structural set of intervals which are a subset of the powerset from 1 to length(x)
#' @param x Data vector: vector of length n > 2
#' @param dec Decay rate: Rate of how fast the different layers will decrease in size. Default is set to sqrt(2).
#' @param minl Minimal theoretical interval-length: The minimal size of intervals to be considered. Default is set to 2.
#' @param thr Threshold: Only intervals with teststatistic greater equal than this threshold will be considered. Should be greater equal to 0. Default is set to 0.
#' @param par Cutoff parameter: Cutoff between full search and optimistic search. Intervals with lower length than the cutoff will be evaulated by full search, the rest will be evaluated by optimistic search, which performs faster for large intervals. Default is set to 40.
#' @param statsonly Binary paramter: If TRUE the function will only calculate and output the teststatics for each interval and omit the changepoint detection. Default is set to FALSE.
#' @param penalty Penalty Type: Two kind of penalty types are possible. If penalty="BIC" the Bayesian information criterion penalty type is used. If penalty="MBIC" a modified BIC criterion is used. Default is set to "BIC".
#' @return
#' A list with the following (except if statsonly=TRUE):
#' \item{Cpts}{Matrix with 6 columns: changepoint found, start of the location where the changpoint was found, end of the location where the changepoint was found, teststatistic to the interval from where the changepoint was taken from, residual sum of squares including the new changepoint and the corresponding criterion penalty}
#' \item{Intervals}{Matrix with 4 columns: start point of the interval, end point of the interval, optimal changepoint for that interval (the optimal splitting point is between output and output+1) and corresponding teststatistic}
#' \item{RSS}{Total residual sum of squares}
#' \item{OptByCrit}{The optimal value for the criterion penalty. If the values fo OptByCrit and RSS coincide then no changepoint is relevant (will not show up if no changepoints are over the threshold)}
#' \item{OptCpts}{All changpoints up to the optimal value for the criterion penalty (will not show up if no changepoints are over the threshold)}
#' @examples
#' \donttest{
#' x=rnorm(10)
#' findcpts(x)
#' }
#' @export
#' @useDynLib ChangePoints


findcpts <-function(x, dec=sqrt(2), minl=2L, thr=mad(diff(x)/sqrt(2)) * sqrt(2 * log(length(x))), par=40L, statsonly=FALSE, penalty="BIC"){
  if (!is.double(x)) {storage.mode(n) <- 'double'}
  if(length(x) < 3){stop("length of x should be at least 3")}
  if (!is.double(dec)) {storage.mode(dec) <- 'double'}
  if (!is.integer(minl)) {storage.mode(minl) <- 'integer'}
  if (!is.double(thr)) {storage.mode(thr) <- 'double'}
  if (!is.integer(par)) {storage.mode(par) <- 'integer'}
  if(statsonly==TRUE){stats=1L}
  else if(statsonly==FALSE){stats=2L}
  else{ stop("statsonly should be TRUE or FALSE")}
  if(penalty=="BIC"){pen=1L}
  else if(penalty=="MBIC"){pen=2L}
  else{ stop("Please input a valid penalty criterion: BIC or MBIC") }
  n <- as.integer(length(x))
  dep <- floor(log(n)/log(dec))
  ilen <- rep(0,dep-1)
  nint <- rep(0L,dep-1)
  for(i in 1:(dep-1)){
    ilen[i] <- n*((1/dec)**i)
    if(ilen[i] < (minl-1)){ dep <- dep-1 }
    else{ nint[i] <- ceiling(n/ilen[i]-1e-8)*2-1 }
  }
  dep <- as.integer(dep)
  nsum <- as.integer(1+sum(nint))
  nint <- as.integer(nint[1:dep])
  ilen <- as.double(ilen[1:dep])
  obj <- .Call("findcptsC",x,n,dec,dep,ilen,nint,nsum,par,stats,pen,thr)
  if(statsonly==TRUE){ return(list(Intervals=cbind(obj$int,obj$value))) }
  else if(obj$t == 0){ return(list(Cpts="none",Intervals=cbind(obj$int,obj$value), RSS=obj$tot)) }
  else if(obj$t2 == 0){ return(list(Cpts=cbind(obj$Cpts2[1:obj$t,],obj$Cpts[1:obj$t,]),Intervals=cbind(obj$int,obj$value),RSS=obj$tot,OptByCrit=obj$tot,OptCpts="none" )) }
  else{ return(list(Cpts=cbind(obj$Cpts2[1:obj$t,],obj$Cpts[1:obj$t,]),Intervals=cbind(obj$int,obj$value),RSS=obj$tot,OptByCrit=obj$Cpts[obj$t2,3],OptCpts=obj$Cpts2[1:obj$t2,1] )) }
}
