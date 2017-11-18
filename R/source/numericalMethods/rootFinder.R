
### rootFinder.R ### 

source("./source/BSformulas.R")
source("./source/SV/heston_wrappers.R")

bisection <- function(target, FUN, ...){
  #
  #   The routine implement the bisection root-finding algorithm to find 
  #   the value x such that FUN(x) = target.  
  #   
  # ARGUMENTS: 
  #   * target = the target value 
  #   * FUN = a function of the form FUN(x) returning a signle value 
  # additional arguments:
  #   * low = lower bound. Default is 0 
  #   * high = upper bound. defaulst is 1.5  
  #   * tollerance = required tollerance. Default is 1.5e-15  
  # VALUE: 
  #   * x 
  
  if(is.na(target)) return(NA)
  
  # additional arguments
  arg <- list(...)
  
  if(is.null(arg$low)) {low=0} else {low=arg$low} 
  if(is.null(arg$high)) {high=1.5} else {high=arg$high}
  if(is.null(arg$tolerance)) {tolerance=1.5e-15} else {tolerance=arg$tolerance}

  x <- 0.5*(low + high);
  y <- FUN(x);
  
  max_loops <- round(log((high-low)/.Machine$double.eps, base=2))
  
  for(i in 1:max_loops){ 
  
    if(y<target) 
      low = x
    if (y>target) 
      high = x
    
    x = 0.5*(low + high)
    y = FUN(x)
    
    if(!(abs(y-target)>tolerance)) {break}
  } 
  
  return(x)

}

BScall_iv <- function(p, S, K, r, q, t, method='bisection', ...){

  FUN <- function(x){BScall_price(S=S, K=K, r=r, q=q, x, t=t)} 

  if(method=='bisection'){
    bisection(FUN, target=p, ...)
  } else if(method=='tangent'){
    print('method not available yet')
  }
}

hBScall_iv <- function(par, S, K, mat, r, q){
  
  FUN <- function(x){
    par[5] <- x
    return(hestonPricer(par, S, K, mat, r, q))
    } 
    bisection(FUN, target=p, ...)
}


test_BS_iv <- function(){
  
  S=100
  K=110
  r=0.05
  q=0.02
  exp=2
  
  IV=0.01
  p_theoretical <- BScall_price(S,K,r,q,IV,exp)
  all.equal(IV, BScall_iv(p_theoretical,S,K,r,q,exp))
  
  IV=4  
  p_theoretical <- BScall_price(S,K,r,q,IV,exp)
  all.equal(IV, BScall_iv(p_theoretical,S,K,r,q,exp, higher=5)) 
  
}

