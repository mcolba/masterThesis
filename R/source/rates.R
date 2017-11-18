
###  rates.R  ###

spotTOforward <- function(rates, t_grid, type = c('contTOcont','contTOsimp')){
  #
  #   The function convert spot rates to forward rates using the flat forward 
  #   convention, i.e. linear on the log of discount factors
  #
  
  # check dimension 
  if(t_grid[1]==0) {stop("t_grid[1]=0")}
  if(length(rates)!=length(t_grid)) {stop("wrong dimension")}
  
  n <- length(t_grid)
  dt <- c(t_grid[1], t_grid[2:n]-t_grid[1:(n-1)])
  
  if(type=='contTOcont'){
    
    f <- c(rates[1], (rates[2:n]*t_grid[2:n] - rates[1:(n-1)]*t_grid[1:(n-1)])/dt[-1])
    
  } else if(type=='contTOsimp'){
    
    r.simp <- (exp(rates*t_grid)-1)/t_grid
    f <-  c(r.simp[1], ((1+r.simp[2:n]*t_grid[2:n])/(1+r.simp[1:(n-1)]*t_grid[1:(n-1)])-1)/dt[-1]) 
    
  }
  return(list(f = f, dt = dt))
}

forwardTOspot <- function(rates, t_grid, type = c('simpTOcont','contTOcont')){
  #
  #   The function convert forward rates to spot rates 
  #
  
  # check dimension 
  if(t_grid[1]==0) {stop("t_grid[1]=0")}
  if(length(rates)!=length(t_grid)) {stop("wrong dimension")}
  
  n <- length(t_grid)
  dt <- c(t_grid[1], t_grid[2:n]-t_grid[1:(n-1)])
  
  if(type=='simpTOcont'){
    
    r.simp <- (cumprod(1+rates*dt)-1)/t_grid
    s <- log(r.simp*t_grid+1)/t_grid
    
  } else if(type=='contTOcont'){
    
    s <- cumsum(rates*dt)/t_grid
    
  }
  return(list(s = s, dt = dt))
}

# Unit testing 
spotTOforward_test <- function(){
  rates <-  seq(0.01, 0.05, length.out = 10)
  t_grid <- seq(0.3, 3, length.out = 10)
  
  dt <- spotTOforward(rates,t_grid, type='contTOcont')$dt
  forward.cont <- spotTOforward(rates,t_grid, type='contTOcont')$f
  forward.simp <- spotTOforward(rates,t_grid, type='contTOsimp')$f
  
  a = exp(rates[length(rates)]*t_grid[length(t_grid)])
  b = prod(exp(forward.cont*dt))
  c = prod(forward.simp*dt+1)
  
  stopifnot(all.equal(a,b) && all.equal(a,c))
  
  spot.cont <- forwardTOspot(forward.cont, t_grid, type='contTOcont')$s
  spot.cont2 <- forwardTOspot(forward.simp, t_grid, type='simpTOcont')$s
  
  stopifnot(all.equal(rates,spot.cont) && all.equal(spot.cont,spot.cont2))
  
}
