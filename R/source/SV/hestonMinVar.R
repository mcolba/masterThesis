
### hestonMinVar.R ### 

source("./source/SV/heston_wrappers.R")

hestonVega <- function(par, S, K, mat, r, q, h=NULL){
  # 
  #   The routine approzimate the vega (derivatives w.r.t. v_0) of an European call 
  #   option under the Heston mode via finite difference.  
  #
  # ARGUMENTS: 
  #   * par = [5]-vector of heston parametes in the order: k, v_inf, volofvol
  #           rho, and v_0.     
  #   * S = underlying spot price 
  #   * K = strike 
  #   * mat = expiry 
  #   * r = continously conpounded (spot) rate
  #   * q = continously conpounded (spot) dividend yield 
  #   * h = finite difference increment
  #
  # VALUE: 
  #   * The function returns vega. 
  # 

  v0_ini <- par[5]

  # set h is not specified by the user 
  if(is.null(h)){
    if (-0.0001<v0_ini & v0_ini<0.0001)
      h = 5E-10
    else
      h = abs(2E-16^(1/3)*v0_ini);    
  }
  
  # central finite difference approx. 
  par[5] <- v0_ini + 0.5*h
  f1 <- hestonPricer(par, S, K, mat, r, q)
  par[5] <- v0_ini - 0.5*h
  f0 <- hestonPricer(par, S, K, mat, r, q)

  return((f1-f0)/h)
  
} 

hestonMinVarDelta <- function(par, S, K, mat, r, q){
  # 
  #   The routine compute the minimum variance (MV) delta for an European call 
  #   option under the Heston mode: delta_MV = delta + rho*volOfVol*vega/S
  #
  # ARGUMENTS: 
  #   * par = [5]-vector of heston parametes in the order: k, v_inf, volofvol
  #           rho, and v_0.     
  #   * S = underlying spot price 
  #   * K = strike 
  #   * mat = expiry 
  #   * r = continously conpounded (spot) rates 
  #   * q = continously conpounded (spot) dividend yields 
  #   * h = finite difference increment
  #
  # VALUE: 
  #   * The function returns the MV delta 
  # 
  
  rho <- par[4]
  volvol <- par[3]
  
  delta <- hestonDelta(par, S, K, mat, r, q) 
  vega <- hestonVega(par, S, K, mat, r, q)
  
  minvar <- delta + rho*volvol*vega/S
  
  return(minvar)  
}
