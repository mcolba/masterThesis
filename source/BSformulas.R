
###  BSformulas.R  ###

BScall_price <- function(S, K, r, q, vol, t){
  
  moneyness = log(S/K)
  std = vol*sqrt(t)
  d1 = (moneyness+(r-q+0.5*vol^2)*t)/std
  d2 = d1 - std 
  return(S*exp(-q*t)*pnorm(d1)-K*exp(-r*t)*pnorm(d2))
  
}

BScall_price2D <- function(spot, strikes, rates, dividends, iv, expiries, ...){
  # 
  # 
  
  ## assert... 
  
  moneyness = log(spot/strikes)
  std <- sqrt(expiries) * iv
  v1 <- (rates - dividends + 0.5*iv^2)*expiries
  d1 <- t(moneyness + t(v1))/std
  d2 <- d1 - std 
  v2 <- spot*exp(-dividends*expiries)*pnorm(d1)
  v3 <- exp(-rates*expiries)*pnorm(d2) %*% diag(strikes)
  
  return(v2-v3)
  
}

BScall_delta <- function(S, K, r, q, vol, t){
  
  moneyness = log(S/K)
  std = vol*sqrt(t)
  d1 = (moneyness+(r-q+0.5*vol^2)*t)/std
  return(exp(-q*t)*pnorm(d1))
  
}

BScall_vega <- function(S, K, r, q, vol, t){
  
  moneyness = log(S/K)
  std = vol*sqrt(t)
  d1 = (moneyness+(r-q+0.5*vol^2)*t)/std
  return(S*exp(-q*t) * sqrt(t)*dnorm(d1))

}

BScall_gamma <- function(S, K, r, q, vol, t){
  
  moneyness = log(S/K)
  std = vol*sqrt(t)
  d1 = (moneyness+(r-q+0.5*vol^2)*t)/std
  return(exp(-q*t)*dnorm(d1)/(S*std))
  
}

BScall_speed <- function(S, K, r, q, vol, t){
  
  gamma <- BScall_gamma(S, K, r, q, vol, t)
  moneyness = log(S/K)
  std = vol*sqrt(t)
  d1 = (moneyness+(r-q+0.5*vol^2)*t)/std
  return(-gamma*(d1/std+1)/S)
  
}

callTOput <- function(call, S, K, r, q, t){
  
  Frw <- S*exp((r-q)*t) 
  tmp <- exp(-r*t)*Frw - matrix(exp(-r*t))%*%K
  return(call - tmp)
}

putTOcall <- function(put, S, K, r, q, t){
  
  Frw <- S*exp((r-q)*t) 
  tmp <- exp(-r*t)*Frw - matrix(exp(-r*t))%*%K
  return(put + tmp)
  
}
  
test_putCallParity <- function(){
  
  S = 100 
  K = 130
  r = 0.05
  q = 0.02
  sigma = 0.25
  t = 1
  
  # Theoretical from www.fintools.com : 
  # c = 2.6055	
    # Delta = 0.2064
  # p = 28.2455
	
  # compute cal price
  c <- BScall_price(S,K,r,q,sigma,t)
  all.equal(c, 2.6055, tolerance=1e-4)

  # compute delta 
  delta <- BScall_delta(S,K,r,q,sigma,t)
  all.equal(delta, 0.2064, tolerance=1e-4)
  
  # compute put price w/ put-call parity 
  p <- callTOput(c,S,K,r,q,t)
  all.equal(p, 28.2455, tolerance=1e-4)
  
  # back to call price 
  c2 <- putTOcall(p,S,K,r,q,t)
  all.equal(c,c2)
  
}