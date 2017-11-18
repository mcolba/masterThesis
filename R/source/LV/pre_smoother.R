
###  pre_smoother.R  ###


# load/install packages 
if (!require("fields")) install.packages("fields") # for Thin Plate Spline
# include R files 
source("./source/BSformulas.R")
source("./source/LV/mesh.R")


TPSpreSmoother <- function(iv, expiries, strikes, spot, rates, dividends, t_grid.new = NULL, 
                           k_grid=NULL, forward=NULL,  add_t = NULL, alpha=NULL, 
                           n.points = 35, ...){
  # 
  # The function perform thin plate smoothing as sudgested in Fengles (2009), 
  # p. 422 algo (i). Convert from strike grid to forward moneyness grid. 
  #
  # ARGUMENTS: 
  # * iv = [NxM] matrix of implied volatilities  
  # * expiries = [1xN] vector of time to matutities in yeears 
  # * strkes = [1xM] vector of strikes or moneyness values 
  # * rates = [1xN] vector of contnously conpounded zero rates 
  # * dividends = [1xN] vector of contnously conpounded dividend yields 
  # * spot = index spot value 
  # * k_grid = vector containing the point in time where the greed will be evaluated,
  #            default is k_grid = strikes. 
  # 
  # VALUE: 
  # * the function returns a list containing the following objects: 
  #   * smoothPrices
  #
  
  # compute forward prices if missing 
  if(is.null(forward)){
    forward <- spot*exp((rates-dividends)*expiries)
  }
  # compute forward moneyness (k) 
  tmp <- expand.grid(forward, strikes)
  moneyness <- tmp[,2]/tmp[,1] 
  # default k_grid (30 equally spaced)
  if(is.null(k_grid)){
    if(is.null(alpha)){
      k_grid <- seq(min(moneyness)-0.05, max(moneyness)+0.05, length.out=n.points)
    } else {
      k_grid <- hyperbolicMesh(min(moneyness)-0.05, max(moneyness)+0.05, c=1, 
                               n=n.points, alpha)$x
    }
  }

  # rename expiries 
  t <- expiries
  
  # trasform IVS to Total Variance 
  totalVar <- as.vector(iv^2*t)
  
  xy <- expand.grid(t, strikes) # create d.f. of coordinates
  xy[,2] <- moneyness
  
  # fit the thin plate spline
  tpsFit <- Tps(xy, totalVar, lambda=1e-12, m = 3) 

  #  set the t_grid use in the prediction and interpolate rates  
  if(!is.null(t_grid.new)){
    idx <- order(t_grid.new) 
    t_grid <-  t_grid.new[idx]
    # Row interpolation (aka flat forward)
    dividends <- approx(x = t, y = dividends*t, xout = t_grid)$y/t_grid
    rates <- approx(x = t, y = rates*t, xout = t_grid)$y/t_grid
    forward <- spot*exp((rates-dividends)*t_grid)
  } else {t_grid = t}
  
  # add additional time points
  if(!is.null(add_t)){
    idx <- order(c(t_grid, add_t)) 
    t_grid <-  c(t_grid, add_t)[idx]
    # Row interpolation (aka flat forward)
    dividends <- approx(x = t, y = dividends*t, xout = t_grid)$y/t_grid
    rates <- approx(x = t, y = rates*t, xout = t_grid)$y/t_grid
    forward <- spot*exp((rates-dividends)*t_grid) 
  } 
  
  # new grid 
  new_xy <- expand.grid(expiry = t_grid, k = k_grid)
  # evaluate on the new grid 
  totVar_smooth <- predict(tpsFit, x = new_xy)
  iv_smooth <- sqrt(totVar_smooth/new_xy[,1])
  
  # prices 
  q <- expand.grid(dividends, k_grid)[,1]
  p <- spot * exp(-q * new_xy[,1]) * (BScall_price(S = 1, K = new_xy[,2], r = 0, q = 0, 
                                                   vol = iv_smooth, t = new_xy[,1]))
  # from vector to matrix 
  smoothPrices = matrix(p, nrow = length(t_grid), ncol = length(k_grid))
  smoothIV <- matrix(iv_smooth, nrow = length(t_grid), ncol = length(k_grid))
  
  #forw <- expand.grid(forward, k_grid)[,1]
  #K <- new_xy[,2]*forw
  
  #df <- BScall_price(spot, K, r, q, iv_smooth, new_xy[,1]) - p
  #max(abs(df))

  # output 
  out = list(p = smoothPrices, 
             iv = smoothIV, 
             xy = new_xy, 
             k_grid = k_grid,
             t_grid = t_grid,
             rates = rates,
             dividends = dividends,
             forward = forward)

  return(out)
  
}







