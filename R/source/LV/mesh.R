
###  mesh.R  ###


# load/install packages 
if(!require("fields")) install.packages("fields") # for Tps 
# include R files 
source("./source/numerical/fdm.R") # for diffMat 


hyperbolicMesh <- function(min, max, c, n, alpha=0.1, x_star=NULL){
  #
  # The function create a grid with points centered around c. See White R., 2013, 
  # 'Numerical Solutions to PDEs with Financial Applications' for more details. 
  #
  # ARGUMENTS:
  # * min = value of the first grid point.                   
  # * max = value of the last grid point.
  # * c = center.  
  # * n = length of the grid. 
  # * alpha = parameter to control the density of points around c. A smaller alpha
  #   will cause the points to be more centered around c. Default is 0.1. 
  # * x_star = point to include in the grid (optional). 
  # 
  # VALUE: 
  # * the function returns a list containing the following objects: 
  #   * $x = grid points 
  #   * $idx = index corresponding to the x_star grid point, 
  #     such that x[idx] = x_star. 
  
  n = n-2 # dont count the extreme points
  
  # create n eqaully spaced z points on [0,1]
  z <- (0:(n+1))/(n+1)
  
  beta <- alpha*(max - min)
  delta <- asinh((min-c)/beta) 
  gamma <- asinh((max-c)/beta)-delta
  
  if(!is.null(x_star)){
    y_star <- (asinh((x_star-c)/beta) - delta)/gamma
    idx <- round(y_star*(n+1)) # index in the 0,...,N notation 
    z_star <- idx/(n+1)
    # fit quadratic a curve 
    A <- rbind(c(0,0,1),c(z_star^2,z_star,1),c(1,1,1)) 
    abc <- solve(A)%*%matrix(c(0,y_star,1))
    # predict new z, which include y_star 
    z <- abc[1]*z^2+abc[2]*z+abc[3] 
  }
  
  # find grid points  
  x <- c + beta*sinh(gamma*z+delta)
  
  # check x_star 
  if(!is.null(x_star)){
    if(!all.equal(x[idx+1],x_star)) stop('x_star do not coincide') 
  } else {idx=NULL}
  
  return(list(x=x,idx=idx+1))
  
}

tpsInterpolation <- function(y, newX, x_mesh, x2_mesh = NULL){
  
  if(is.null(x_mesh)){
    x <- x_mesh
  } else {
    x <- expand.grid(x_mesh, x2_mesh)
    newX <- expand.grid(newX, x2_mesh)
  }
  
  tpsFit <- Tps(x, z, lambda = 1e-2)
  out <- predict(tpsFit, x = newX)
  
  return(out)
}

flatDeltaExtrap <- function(iv, k_grid){
  
  # check order 
  idx <- order(k_grid)
  k_grid <- k_grid[idx]
  iv <- iv[,idx]
  
  # dim 
  nk <- length(k_grid)
  nt <- nrow(iv)
  dx.l <- k_grid[2]-k_grid[1]
  dx.u <- k_grid[nk]-k_grid[nk-1]  
  
  # output arrays 
  k.out <- c(k_grid[1]-dx.l, k_grid, k_grid[nk]+dx.u)
  iv.out <- cbind(rep(NA,nt),iv,rep(NA,nt))

  # derivatives 
  iv_k <- t(apply(iv, 1, function(x){diffMat(k_grid, 1)%*%x}))
  iv_kk <- t(apply(iv, 1, function(x){diffMat(k_grid, 2)%*%x}))
  
  # left extrapol 
  iv.out[,1] <- iv.out[,2] - dx.l*iv_k[,1] + dx.l^2*iv_kk[,1]/4
  # right extrapol 
  iv.out[,nk+2] <- iv.out[,nk+1] + dx.u*iv_k[,nk] + dx.u^2*iv_kk[,nk]/4
  
  return(list(iv=iv.out, K=k.out))
}  

shiftIV <- function(data,t,days){
  
  newData <- data[[t+days]]
  newExp <- data[[t]]$expiries-days/250
  
  for(i in 1:ncol(newData$iv)){
    newData$iv[,i] <- approx(x = newData$expiries, 
                             y = newData$iv[,i], 
                             xout = newExp, rule=2)$y
  }
  newData$rates <- approx(x = newData$expiries, 
                          y = newData$rates*newData$expiries, 
                          xout = newExp, rule=2)$y/newExp
  newData$dividends <- approx(x = newData$expiries, 
                          y = newData$dividends*newData$expiries, 
                          xout = newExp, rule=2)$y/newExp
  newData$expiries <- newExp  
    
  return(newData)
}

