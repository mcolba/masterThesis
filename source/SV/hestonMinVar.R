
### hestonMinVar.R ### 

source("./source/SV/heston_wrappers.R")



hestonVega <- function(par, S, K, mat, r, q, h=1e-6){
  
  v0_ini <- par[5]
  
  # using finite difference 
  par[5] <- v0_ini + 0.5*h
  f1 <- hestonPricer(par, S, K, mat, r, q)
  par[5] <- v0_ini - 0.5*h
  f0 <- hestonPricer(par, S, K, mat, r, q)

  return((f1-f0)/h)
  
} 

hestonMinVarDelta <- function(par, S, K, mat, r, q){
  
  rho <- par[4]
  volvol <- par[3]
  
  delta <- hestonDelta(par, S, K, mat, r, q) 
  vega <- hestonVega(par, S, K, mat, r, q)
  
  minvar <- delta + rho*volvol*vega/S
  
  return(minvar)  
}
