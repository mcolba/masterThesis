
#  cubicSpline.R 

# 
#  c.spline_predict  
# 
c.spline_predict <- function(newK, u, g, gamma){
  #
  # Given a set of natural cubic splines, the function predict values, first order derivatives, 
  # and second order derivatives.  
  #
  # ARGUMENTS:
  # * newK = points where to evaluate the spline. Can be either a singol value or a matrix. 
  # * u = [NxM]-matrix of knots. N is the time-dimension and M the K-dimension. i.e. each row 
  #       correspond to a natural cubic spline. 
  # * g = [MxM]-matrix of values, where N is the time-dimension and M the K-dimension. 
  # * gamma = [NxM]-matrix of second order derivatives, where N is the time-dimension 
  #           and M the K-dimension.
  # 
  # VALUE: 
  # * the function returns a list containing the following objects: 
  #   * $g = values  
  #   * $delta = first order derivatives 
  #   * $gamma = second order derivatives 
  # 
  
  # check arguments dim. 
  if(is.vector(g)){g = matrix(g, nrow=1)}
  if(is.vector(u)){u = matrix(u, nrow=1)}
  if(is.vector(gamma)){gamma = matrix(gamma, nrow=1)}
  
  if(!identical(dim(u), dim(g), dim(gamma))){
    stop("wrong dimension")
  }
  
  nt <- nrow(g)
  nk <- ncol(g)
  
  if(is.vector(newK)){
    m <- length(newK)
    newK <- matrix(rep(newK, nt), nt, m, byrow = T)
  } 
  
  # from value-second derivative to a, b, c, d representation
  h <- matrix(u[,2:nk]-u[,1:(nk-1)], nt, nk-1)
  m <- nk+1
  # inizialize  
  a <- matrix(NA, nt, m)
  b <- matrix(NA, nt, m)
  c <- matrix(NA, nt, m)
  d <- matrix(NA, nt, m)
  # compute constants, see Fangler (2006)
  a[,1] <- g[,1]
  a[,2:m] <- g
  b[,2:(m-1)] <- (g[,2:nk]-g[,1:(nk-1)])/h - h*(2*gamma[,1:(nk-1)]+gamma[,2:nk])/6
  b[,1] <- b[,2]
  b[,m] <- (g[,nk]-g[,(nk-1)])/h[,(nk-1)] + h[,(nk-1)] * gamma[,nk-1] / 6 
  c[,2:(m-1)] <- 0.5 * gamma[,1:(nk-1)]
  c[,c(1,m)] <- 0
  d[,2:(m-1)] <- (gamma[,2:nk] - gamma[,1:(nk-1)])/(6*h)
  d[,c(1,m)] <- 0   
  
  # initialize output  
  g.o <- matrix(NA, nt, ncol(newK))
  delta.o <- matrix(NA, nt, ncol(newK)) 
  gamma.o <- matrix(NA, nt, ncol(newK)) 
  u.o <- matrix(NA, nt, ncol(newK))
  
  # loop 
  for(t in 1:nt){
    for(i in 1:length(newK[t,])){
      for(j in 1:(nk-1)){
        
        if ((u[t,j]<=newK[t,i]) && (newK[t,i]<u[t,j+1])){
          # intrapolation 
          du <- newK[t,i]-u[t,j]
          g.o[t,i] <- d[t,(j+1)]*du^3 + c[t,(j+1)]*du^2 + b[t,(j+1)]*du +a[t,(j+1)] 
          delta.o[t,i] <- 3*d[t,(j+1)]*du^2 + 2*c[t,(j+1)]*du + b[t,(j+1)] 
          gamma.o[t,i] <- 6*d[t,(j+1)]*du + 2*c[t,(j+1)]
          break
          
        } else if(u[t,nk] <= newK[t,i]){
          # extrapolation 
          du <- newK[t,i] - u[t,nk]
          g.o[t,i] <- b[t,m]*du + a[t,m] 
          delta.o[t,i] <- 3*d[t,m]*du^2 + 2*c[t,m]*du + b[t,m] 
          gamma.o[t,i] <- 6*d[t,m]*du + 2*c[t,m] 
          warning(paste("values at grid point i =", t, " j =", i, "was extrapolated"))
          break 
          
        } else if(newK[t,i] < u[t,1]){
          # extrapolation 
          du <- newK[t,i] - u[t,1]
          g.o[t,i] <- b[t,1]*du +a[t,1] 
          delta.o[t,i] <- 3*d[t,1]*du^2 + 2*c[t,1]*du + b[t,1] 
          gamma.o[t,i] <- 6*d[t,1]*du + 2*c[t,1] 
          warning(paste("values at grid point i =", t, " j =", i, "was extrapolated"))
          break
        }
      }
    }    
  }
  return(list(g = g.o, delta = delta.o, gamma = gamma.o))
} 

#
# evaluate 
#
evalSplin <- function(point, u, g, gamma){
  #
  # 
  # 
  
  if(any(is.na(g)) || any(is.na(gamma))){
    stop("NA values")
  }
  
  n <- length(u)
  if (point <= min(u)){
    dg = (g[2]-g[1])/(u[2]-u[1]) - 1/6*(u[2]-u[1])*gamma[2]
    y = g[1] - (u[1] - point)*dg
    return(y)
  } else if (point>max(u)){
    dg = (g[n]-g[n-1])/(u[n]-u[n-1]) + 1/6*(u[n]-u[n-1])*gamma[n-1]
    y = g[n] + (point - u[n])*dg
    return(y)
  } else {
    for (i in 1:(n-1)){
      if((point > u[i]) && (point <= u[i+1])){
        h = u[i+1] - u[i]
        y = ((point-u[i])*g[i+1] + (u[i+1]-point)*g[i])/h - 1/6*(point-u[i])*(u[i+1]-point) * 
          ( (1+(point-u[i])/h)*gamma[i+1] + (1+(u[i+1]-point)/h)*gamma[i] )
        return(y)
      }
    }    
  }
}


c.spline_delta <- function(i, u, g, gamma){

  n = length(u)
  if(i==0) {i = 1}   
  
  if((1<=i) && (i<n)) {
    
    b = (g[i+1]-g[i])/(u[i+1]-u[i]) - 1/6*(u[i+1]-u[i])*(2*gamma[i] + gamma[i+1])
    return(b)
    
  } else if(i == n){
    
    b = (g[n]-g[n-1])/(u[n]-u[n-1]) + 1/6*(u[n]-u[n-1])*gamma[n-1] # gamma(n-1) or (n-2) ??? 
    return(b)
    
  } else (stop("index out of bounds"))
  
}

#
#
getDeltaGrid <- function(u, g, gamma){
  
  nt <- nrow(u)
  nk <- ncol(u)
  out <- matrix(NA, nt, nk)
  
  for (i in 1:nt){
    for (j in 1:nk){
      out[i,j] <- c.spline_delta(j, u[i,], g[i,],gamma[i,])
    }
  }
  return(out)
} 


# test 
#g = cs$p
#u = cs$K
#gamma = cs$gamma

#new <- c.spline_predict(u, u, g, gamma)
#quantile(new$g - g)
#quantile(new$gamma - gamma)
#quantile(new$delta - getDeltaGrid(u, g, gamma))

#c.spline_predict(1250, u[1,], g[1,], gamma[1,])$g - evalSplin(1250, u[1,], g[1,], gamma[1,])
#c.spline_predict(1000, u[1,], g[1,], gamma[1,])$g - evalSplin(1000, u[1,], g[1,], gamma[1,])
#c.spline_predict(2000, u[1,], g[1,], gamma[1,])$g - evalSplin(2000, u[1,], g[1,], gamma[1,])



