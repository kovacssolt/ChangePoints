#' Changpoint Finder
#'
#' Function that finds change points in given data set x by evaluating a set of seeded intervals and performing GREEDY SELECTION.
#' @param x Data: numeric vector of length n > 2.
#' @param dec Decay rate: Rate of how fast the different layers will decrease in size. Default is set to sqrt(2). In general, this should be a numeric value in (1,2].
#' @param minl Minimal theoretical interval length in seeded intervals: The minimal size of intervals to be considered. Default is set to 2. In general, this should be an integer in [2,length(x)].
#' @param thr Threshold: Only intervals with test statistic (CUSUM) greater or equal than this threshold will be considered. Default is set based on asymptotics. In general, this should be a numeric value in [0, infinity).
#' @param par Cutoff parameter: Cutoff between full grid search and (naive) optimistic search. Intervals with length lower than the cutoff will be evaulated by full grid search, the rest will be evaluated by (naive) optimistic search. Optimistic search is faster for large intervals (roughly from with length bigger than 40), but does not guarantee to find the best split point. Default is set to length(x) which means that full grid search is performed in all intervals.
#' @param statsonly T/F paramter: If TRUE the function will only calculate and output the test statics (CUSUM) for each interval and omit the selection. Default is set to FALSE, and in that case greedy selection is performed for all intervals with a test statistics value above thr. 
#' @param penalty Penalty Type: Two kind of penalty types are possible. If penalty="BIC" the Bayesian information criterion penalty type is used. If penalty="MBIC", the modified BIC criterion is used. Default is set to "BIC". 
#' @return
#' A list with the following (except if statsonly=TRUE):
#' \item{Res}{Matrix with 6 columns: found change point; start of interval in which the change point was found; end of the interval in which the change point was found; test statistic (CUSUM) of the (seeded) interval from where the changepoint was taken from (hence, this is a value that can be found in "Intervals" as well); overall sum of squared errors when including the new changepoint additionally to previous ones; and the value of the information criterion corresponding to the segmentation.}
#' \item{Intervals}{Matrix with 4 columns: start point of the (seeded) interval; end point of the (seeded) interval; found change point withing the interval given by the first two entries; and corresponding test statistic (CUSUM)}
#' \item{RSS}{Total residual sum of squares of the model without any change points.}
#' \item{OptInfCrit}{The optimal (minimal) value for the model selection criterion along the solution path. Note: If none of the intervals had a test statistic above the specified threshold (such that the Res matrix has zero rows) or the model with no change point has a smaller corresponding information criterion than models with change points, then this value will coresspond to the information criterion of the model with no change points.)}
#' \item{OptCpts}{Change points from the optimal segmentation from the solution path (according to the chosen information criterion). Note that the solution path is only calculated for intervals/candidates with a test statistic (CUSUM) above the chosen threshold.}
#' @examples
#' \donttest{
#' x=rnorm(10)
#' findcpts(x)
#' }
#' @export
#' @useDynLib ChangePoints


findcpts <-function(x, dec=sqrt(2), minl=2L, thr=1.3/2*mad(diff(x)/sqrt(2)) * sqrt(2 * log(length(x))), par=length(x), statsonly=FALSE, penalty="BIC"){
  if(!is.numeric(x) | !is.vector(x)) {stop("x is not a numeric vector")}
  if(length(x) < 3){stop("length of x should be at least 3")}
  if(!is.double(x)) {storage.mode(x) <- 'double'}
  
  if(!is.numeric(dec) | !is.vector(dec)) {stop("dec is not a numeric value")}
  if(length(dec)!=1 | any(dec) <= 1 | any(dec) > 2){stop("dec should be a single numeric value in (1,2]")}
  if(!is.double(dec)) {storage.mode(dec) <- 'double'}

  if(!is.numeric(minl) | !is.vector(minl)) {stop("minl is not an integer value")}
  if(length(minl)!=1 | any(minl) <= 1 | any(minl) > length(x)){stop("minl should be a single integer value bigger (or equal) than 2 and less or equal than the length of x")}
  if(!is.integer(minl)) {storage.mode(minl) <- 'integer'}

  if(!is.numeric(thr) | !is.vector(thr)) {stop("thr is not a numeric value")}
  if(length(thr)!=1 | any(thr) < 0){stop("thr should be a single non-negative value")}
  if (!is.double(thr)) {storage.mode(thr) <- 'double'}
  
  if(!is.numeric(par) | !is.vector(par)) {stop("par is not an integer value")}
  if(length(par)!=1 | any(par) < 0 | any(par) > length(x)){stop("par should be a single integer value between 1 and length of x")}
  if (!is.integer(par)) {storage.mode(par) <- 'integer'}
  
  if(statsonly==TRUE){stats=1L}
  else if(statsonly==FALSE){stats=2L}
  else{stop("statsonly should be TRUE or FALSE")}
    
  if(penalty=="BIC"){pen=1L}
  else if(penalty=="MBIC"){pen=2L}
  else{stop("Please input a valid penalty criterion: BIC or MBIC")}
    
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
  colnames(obj$Cpts2) <- c("cpt", "st", "end")
  colnames(obj$Cpts)  <- c("CUSUM_original", "SSE", "InfCrit")
  colnames(obj$int)   <- c("st", "end", "cpt")
  colnames(obj$val)   <- c("CUSUM") 
  
  if(statsonly==TRUE){ 
    return(list(Intervals=cbind(obj$int,obj$value))) 
  }
  else if(obj$t == 0){
    warning("There was no interval with a test statistic above chosen threshold. Check if thr was set properly.") 
    return(list(Res=matrix(nrow=0,ncol=6), Intervals=cbind(obj$int,obj$value), RSS=obj$tot, OptInfCrit=n/2*log(obj$tot/n), OptCpts=integer()))
  }
  else if(obj$t2 == 0){
    return(list(Res=cbind(obj$Cpts2[1:obj$t,], obj$Cpts[1:obj$t,]), Intervals=cbind(obj$int,obj$value), RSS=obj$tot, OptInfCrit=n/2*log(obj$tot/n), OptCpts=integer()))
  }
  else{
    return(list(Res=cbind(obj$Cpts2[1:obj$t,], obj$Cpts[1:obj$t,]), Intervals=cbind(obj$int,obj$value), RSS=obj$tot, OptInfCrit=obj$Cpts[obj$t2,3], OptCpts=obj$Cpts2[1:obj$t2,1]))
  }
}
