
###  extrapolation.R  ###


# include R files 
source("./source/numerical/fdm.R") # for diffMat 


flateningEtrapol <- function(x, grid, internal, method='exponential', lambda=NULL){
  #
  #
  # DETAILS: 
  # * Method 'exponential' uses a convex-decreasing first-order derivative, the lambda 
  #   parameter controls the spead at which delta converges to zero. The higher the 
  #   lambda the faster is the convergence, a value of 4.6 will produce adelta at the end 
  #   of the new grid equal to 1% of the boundry delta from the initial grid. lambda=0 
  #   produces a constant delta.
  # * Method 'power' uses a monotone-decreasing first-order derivative, going from the 
  #   the boundry delta of  the initial grid to zero. The lambda parameter controls 
  #   the curvature: 
  #   * lambda > 1 ->  concave
  #   * lambda = 1 -> linear 
  #   * lambda > 1 -> convex 
  
  # dim 
  N = length(grid)
  n = sum(internal)
  
  # checks 
  if(n!=ncol(x)){stop('wrong dimension')}
  if(!all(order(grid)==(1:N))){stop('wrong order')}
  
  # set lambda if null 
  if(is.null(lambda)){
    if(method=='exponential'){lambda=0} # constant 
    if(method=='power'){lambda=1} # linearly decreasing 
  }
  
  # output matrix 
  x.out <- matrix(NA,nrow(x),N)
  
  # indexing 
  external.low <- grid[internal][1]>grid
  external.upp <- grid[internal][n]<grid
  
  # compute 1st order derivatives 
  D1 <- diffMat(grid[internal])$mat
  deltas <- D1%*%t(x)
  delta.l <- deltas[1,]
  delta.u <- deltas[n,]
  
  # extrap. points dim 
  nu <- sum(external.upp)
  nl <- sum(external.low)
  u.bound <- nl + n
  l.bound <- nl + 1
  
  # dx 
  
  dx.u <- grid[external.upp] - c(grid[u.bound],grid[external.upp][-nu])
  cum_dx.u <- grid[external.upp] - grid[u.bound]
  if(nl>1){
    dx.l <- c(grid[l.bound], grid[external.low][nl:2]) - grid[external.low][nl:1]
    cum_dx.l <- grid[l.bound]  - grid[external.low][nl:1] 
  } else (dx.l <- cum_dx.l <- grid[l.bound]  - grid[external.low][nl:1])
  
  if(method=='exponential'){
    #method 1
    discount.l <- 1/exp((1:nl)/nl*lambda)
    discount.u <- 1/exp((1:nu)/nu*lambda)
  } else if(method=='power'){
    # method 2 
    discount.l <- (1-((1:nl)/nl)^lambda)
    discount.u <- (1-((1:nu)/nu)^lambda)  
  } else{stop('unknown method')}
  
  # lower extrapolation 
  deltas <- matrix(delta.l)%*%discount.l
  Ddx <- t(apply((t(deltas)*dx.l),2,cumsum))
  x.out[,external.low] <- x[,1] - Ddx[,nl:1]
  
  # upper extrapolation 
  deltas <- matrix(delta.u)%*%discount.u
  Ddx <- t(apply((t(deltas)*dx.u),2,cumsum))
  x.out[,external.upp] <- x[,n] + Ddx
  
  # add internal points 
  x.out[,internal] <- x 
  
  return(x.out)
}

benaim.fit <- function(mu, K, p, p_k, p_kk, type='call'){
  #
  # Fit the function proposed by Benaim et al. to the cut-off call 
  # (or put) option value and its first and second order derivatives.    
  #
  
  if(type=='call'){
    
    # Solve Hx = g, for x = [b,c]'  
    # 2nd equation terms 
    h1 = c(1/(K^2), 2/(K^3))
    g1 = -p_k/p - mu/K 
    # 3rd equation terms 
    h2 = c(-p_k/(K^2) + 2*p/(K^3), 6*p/(K^4) - 2*p_k/(K^3)) 
    g2 = p_k*mu/K - mu*p/(K^2) + p_kk 
    # write the system in matrix form 
    H = rbind(h1, h2)
    g = matrix(c(g1,g2))
    # solve 
    sol = solve(H)%*%g 
    
    c = sol[2]
    b = sol[1]
    a = log(p*K^mu) - b/K - c/(K^2)
    
    return(c(a,b,c,mu))
    
  } else if(type=='put'){
    
    # Solve Hx = g, for x = [b,c]'  
    # 2nd equation terms 
    h1 = c(1, 2*K)
    g1 = p_k/p - mu/K 
    # 3rd equation terms 
    h2 = c(p_k, p_k*2*K + 2*p) 
    g2 = p_kk - p_k*mu/K +  p*mu/K^2 
    # write the system in matrix form 
    H = rbind(h1, h2)
    g = matrix(c(g1,g2))
    # solve 
    sol = solve(H)%*%g 
    
    c = sol[2]
    b = sol[1]
    a = log(p/K^mu) - b*K - c*K^2
    
    return(c(a,b,c,mu))
  }
  

}

benaim.predict <- function(K, par, type='call'){
  #
  # predict the tail values for a value or vector of strikes K.  
  #
  
  if(type=='put'){K=1/K}

  a <- par[1]
  b <- par[2]
  c <- par[3]
  mu <- par[4]
  
  return(K^-mu * exp(a+ b/K + c/(K^2)))
  
}


benaimEtrapol <- function(p.call , new.grid, internal, spot, r, q, t_grid,
                          mu.call=3, mu.put=3, c_k.analytic=NULL, c_kk.analytic=NULL){
  
  # x = [M,N]-matrix. Rows correspond to expirations and columns to strikes. 
  # new.grid = k_grid, must have an increasing order. 
  
  if(is.vector(p.call)){ p.call = matrix(p.call,nrow = 1) } 
  if(is.vector(c_k.analytic)){ c_k.analytic = matrix(c_k.analytic) }
  if(is.vector(c_kk.analytic)){ c_kk.analytic = matrix(c_kk.analytic) }
  
  # dim 
  N = length(new.grid)
  n = sum(internal)
  nt = length(t_grid)
  
  # checks 
  if(n!=ncol(p.call)){stop('wrong dimension')}
  if(nt!=nrow(p.call)){stop('wrong dim')}
  if(!all(order(new.grid)==(1:N))){stop('wrong order')}
  
  # output matrix 
  p.out <- matrix(NA,nt,N)
  
  # indexing 
  external.low <- new.grid[internal][1]>new.grid
  external.upp <- new.grid[internal][n]<new.grid
  
  # diff. marices 
  D1 <- diffMat(new.grid[internal],1,accurate.boundry = F)$mat 
  D2 <- diffMat(new.grid[internal],2,accurate.boundry = F)$mat
  
  # cal prices 
  call <- p.call[,n]
  # put prices via put-call parity 
  p.put <- callTOput(p.call, spot, K=new.grid[internal], r, q, t_grid)
  put <- p.put[,1]

  # derivatives 
  if(!is.null(c_k.analytic)){
    # use analytical first order der. 
    call_k <- c_k.analytic[n,]
  } else { call_k <- (D1%*%t(p.call))[n,] } 
  
  if(!is.null(c_kk.analytic)){
    # use analytical second order der.
    call_kk <- c_kk.analytic[n,]
  } else { call_kk <- (D2%*%t(p.call))[n,] }
  # put derivatives 
  put_k <- (D1%*%t(p.put))[1,]
  put_kk <- (D2%*%t(p.put))[1,]
  
  # mat to store tmp extrap put prices 
  put.exterap <- matrix(NA,nt,sum(external.low))
  
  # loop 
  for(i in 1:nt){
    
    # put extrapol
    if(put[i]!=0){
      
      fit <-  benaim.fit(mu=mu.put, K=new.grid[internal][1], p=put[i], p_k=put_k[i], 
                         p_kk=put_kk[i], type='put')
      put.exterap[i,] <- benaim.predict(K=new.grid[external.low], par=fit, type='put')
      
    } else {put.exterap[i,external.low] <- 0}
    
    # call extrapol
    if(call[i]>1e-15)
    {
      # check delta and gamma 
      if(call_k[i]>-1e-15){
        p.out[i,external.upp] <- call[i]
        warning(paste('Negative c_k at time ', i, ' was set to zero'))
        next
      }
      if(call_kk[i]<1e-15){
        call_kk[i]=0
        warning(paste('Negative c_kk at time ', i, ' was set to zero'))
      }
      # fit the extrap. function 
      fit <-  try(benaim.fit(mu=mu.put, K=new.grid[internal][n], p=call[i], p_k=call_k[i], 
                         p_kk=call_kk[i], type='call'), silent=TRUE)
      if(is(fit, "try-error")){
        warning('Benaim.fit did not find a feasable solution (singular matrix) at time ',
                i,'. A flat extrapolation was used instead' )
        p.out[i,external.upp] <- call[i]
        next
      } else {
        pred <- benaim.predict(K=new.grid[external.upp], par=fit, type='call')
        p.out[i,external.upp] <- pred  
        
        # chek that the result is decreasing and finite 
        if(pred[1]<pred[sum(external.upp)]){ 
          warning(paste('Extrapolation failed, increasing pricing at time ', i, 
                        '. A flat extrapolation was used instead'))
          p.out[i,external.upp] <- call[i]
          next
        } 
        is_infORneg <- is.infinite(p.out[i,external.upp]) | p.out[i,external.upp]<0
        p.out[i,external.upp][is_infORneg] <- 0
      }
    } else {p.out[i,external.upp] <- 0}
    
  }
  
  # back to call prices 
  p.out[,external.low] <- putTOcall(put.exterap,spot,new.grid[external.low],r,q,t_grid)
  p.out[,internal] <- p.call
  
  return(p.out)
}


test_benaim <- function(){
  
  # call option cut-off parameters   
  K=140     # strike
  S=100     # spot 
  p=0.22       # call price  
  p_k=-0.035   # delta 
  p_kk=0.005 # gamma
  # fit and predict 
  fit <- benaim.fit(5,K,p,p_k,p_kk)
  tail.K <- seq(140,200,length.out=20) 
  tail.p <- benaim.predict(tail.K, fit)
  # asert 
  all.equal(tail.p[length(tail.p)], 0, tolerance=1e-12)
  all.equal(tail.p[1], p, tolerance=1e-12)
  all(diffMat(tail.K,1)$mat%*%tail.p<=0)
  all(diffMat(tail.K,2)$mat%*%tail.p>=0)
  
}
