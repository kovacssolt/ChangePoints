#' Changepoint Finder with narrow interval prioritization
#'
#' Function that finds changepoints in given data set x over a specific threshold and prioritizing narrow intervals
#'
#' @param x Data vector: vector of length n > 2
#' @param dec Decay rate: Rate of how fast the different layers will decrease in size. Default is set to sqrt(2).
#' @param minl Minimal theoretical interval-length: The minimal size of intervals to be considered. Default is set to 2.
#' @param thr Positive Threshold: Only intervals with teststatistic greater equal than this threshold will be considered. Default is set to 1e-4.
#' @param par Cutoff parameter: Cutoff between full search and optimistic search. Intervals with lower length than the cutoff will be evaulated by full search, the rest will be evaluated by optimistic search, which performs faster for large intervals. Default is set to 40.
#' @return
#' A list with the following:
#' \item{Cpts}{Sorted changepoints found}
#' \item{Intervals}{Matrix with 4 columns: start point of the interval, end point of the interval, optimal changepoint for that interval (the optimal splitting point is between output and output+1) and corresponding teststatistic}
#' @examples
#' \donttest{
#' x=rnorm(10)
#' narrowcpts(x,thr=0.1)
#' }
#' @useDynLib ChangePoints


narrowcpts <-function(x, dec=sqrt(2), minl=2L, thr=1e-6, par=40L){
  if (!is.double(x)) {storage.mode(n) <- 'double'}
  if(length(x) < 3){stop("length of x should be at least 3")}
  if (!is.double(dec)) {storage.mode(dec) <- 'double'}
  if (!is.integer(minl)) {storage.mode(minl) <- 'integer'}
  if (!is.double(thr)) {storage.mode(thr) <- 'double'}
  if(thr <= 0){stop("thr should be positive")}
  if (!is.integer(par)) {storage.mode(par) <- 'integer'}
  n <- as.integer(length(x))
  dep <- floor(log(n)/log(dec))
  ilen <- rep(0,dep-1)
  nint <- rep(0L,dep-1)
  for(i in 1:(dep-1)){
    ilen[i] <- n*((1/dec)**i)
    if(ilen[i] < (minl-1)){ dep <- dep-1 }
    else{ nint[i] <- ceiling(round(n/ilen[i], 14))*2-1 }
  }
  dep <- as.integer(dep)
  nsum <- as.integer(1+sum(nint))
  nint <- as.integer(nint[1:dep])
  ilen <- as.double(ilen[1:dep])
  obj <- .Call("narrowcptsC",x,n,dec,dep,ilen,nint,nsum,par,thr)
  if(obj$t == 0){ return(list(Cpts="none",Intervals=cbind(obj$int,obj$value))) }
  else{ return(list(Cpts=sort(obj$Cpts2[1:obj$t]),Intervals=cbind(obj$int,obj$value))) }
}
