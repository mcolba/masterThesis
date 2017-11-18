
###  Fengler09.R  ### 


# load/install packages 
if (!require("quadprog")) install.packages("quadprog") # for solve.QP 
# include R files 
source("./source/LV/mesh.R")
source("./source/BSformulas.R")


solveQuadprog <- function(p, k_grid, t_grid, rates, dividends, spot, 
                          forward, smooth = 1e-2, sig.space=NULL, ...){
  # 
  # The function solve the quadratic problem defined in Fengles (2009), p. 422 algo (ii):
  #   min_x   -y'x+0.5x'Bx
  #   s.t.    A'x = 0 
  #           + inequality constraints 
  #
  # ARGUMENTS:
  # * p = [N x M] matrix of call prices 
  # * 
  # 
  # VALUE: 
  # * the function returns a list containing the following objects: 
  #   * $g 
  #   * $gamma 
  # 
  
  if(is.null(sig.space)){
    sig.space <- matrix(TRUE,length(t_grid),length(k_grid))
  } 
  
  # dim 
  nt <- length(t_grid) 
  dt <- c(t_grid[1], t_grid[2:nt] - t_grid[1:(nt-1)])
  nk <- length(k_grid)
  
  # check dim 
  if(nk!=ncol(p) | nk!=ncol(sig.space)){ stop('wrong dim') }
  if(nt!=nrow(p) | nt!=nrow(sig.space)){ stop('wrong dim') } 
  
  # impose right order 
  idx <- order(k_grid)
  k_grid <- k_grid[idx]
  p <- p[,idx]; sig.space <- sig.space[,idx]
  
  # compute one step forward rates and forward dividends 
  if (length(rates)==length(dividends) && length(dividends)==nt){
    forwardDividends <- spotTOforward(dividends, t_grid, 'contTOcont')$f
  } else {stop("non conformable arrays")}
  
  u <- matrix(forward) %*% t(matrix(k_grid)) # [nt x nk]-matrix of knots 
  
  g <- matrix(NA, nt, nk) # value 
  gamma <- matrix(NA, nt, nk) # value-second derivative 

  # loop 
  for (j in nt:1){
    
    lambda = smooth
    
    # reduce grid dimension  
    nk <- sum(sig.space[j,])
    u_j <- u[j,sig.space[j,]]
    p_j <- p[j,sig.space[j,]]

    # h <- u[j,2:nk] - u[j,1:(nk-1)] # [1,nk-1]-vector of knots first differences 
    h <- u_j[2:nk] - u_j[1:(nk-1)] # [1,nk-1]-vector of knots first differences 
    
    Q <- matrix(0, nk, nk-2)
    for (i in 2:(nk-1)){
      Q[i-1,i-1] <- 1/h[i-1]
      Q[i, i-1] <- -1/h[i-1]-1/h[i]
      Q[i+1, i-1] <- 1/h[i]
    }    
    
    R <- matrix(0, nk-2, nk-2)
    for (i in 2:(nk-1)){
      R[i-1,i-1] = 1/3*(h[i-1]+h[i]);
      if (i<nk-1){
        R[i-1,i] = 1/6*h[i];
        R[i,i-1] = 1/6*h[i];        
      }
    }
    
    weights <- rep(1,nk) # no weights, as in Fengler (2009)
    # iv <- seq(0.2,0.1,length.out=nk)
    # weights <- abs(1/BScall_vega(S0,u_j,0,0,iv,t_grid[j]))
    
    # y <- matrix(c(weights*p[j,],rep(0,nk-2))) # [2nk-2, 1]
    y <- matrix(c(weights*p_j,rep(0,nk-2))) # [2nk-2, 1]
    
    A <- rbind(Q, t(-R)) 
    B <- rbind( cbind(diag(weights), matrix(0, nk, nk-2)),
                cbind(matrix(0, nk-2, nk), lambda*R) )
    
    # CONSTRAINTS, see algo. (ii) 
    # x = (g_1,...,g_n,gamma_2,...,gamma_n-1)' is a [2nk-2, 1]-vector and we need to 
    # create the matrices LHS [2k-2, ..] and rhs [.., 1] such that C'x >= d.       
    
    # (*)
    if(j==nt){
      nocs_lhs <- c(-1, rep(0, 2*nk-3)) # general no-arbitrage (eq. 6 or 18)
      nocs_rhs <- -spot*exp(-dividends[j]*t_grid[j]) # [1]
    } else {
      nocs_lhs <- cbind(diag(-1,nk,nk), matrix(0,nk,nk-2)) # no-calendar spread (prop. 2.1)
      # nocs_rhs <- -g[j+1,] * exp(forwardDividends[j+1]*dt[j+1]) # [nk, 1] 
      nocs_rhs <- -g[j+1,sig.space[j,]] * exp(forwardDividends[j+1]*dt[j+1]) # [nk, 1] 
    }
    
    # LHS'
    LHS_t <- rbind(
      t(A), # EQUALITY for getting a cubic natural spline 
      cbind(matrix(0, nk-2, nk), diag(rep(1,nk-2))), # convex in K (eq. 16)
      # c(-1/h[1], 1/h[1], rep(0,nk-2), -h[1]/6, rep(0,nk-3)), # non-increasing in K (eq. 4 & 17)  
      # or: 
      c(-1/h[1], 1/h[1], rep(0,nk-2), rep(0,nk-2)),
      # c(rep(0,nk-2), 1/h[nk-1], -1/h[nk-1], rep(0,nk-3), -h[nk-1]/6), # non-increasing in K (eq. 4 & 17)
      # or: 
      c(rep(0,nk-2), 1, -1, rep(0,nk-2)),
      nocs_lhs, # (*)
      c(1, rep(0, 2*nk-3)), # general no-arbitrage (eq. 6 & 18)
      c(rep(0, nk-1), 1, rep(0, nk-2)) # g_n non-negative (e.q 6 & 18)   
    )
    # RHS 
    rhs <- rbind(
      matrix(0, nk-2, 1), # [nk-2,1]
      matrix(0, nk-2, 1), # [nk-2,1]
      matrix(-exp(-rates[j]*t_grid[j])), # [1,1]
      matrix(0),# [1,1]
      matrix(nocs_rhs), # (*)
      # matrix(spot*exp(-dividends[j]*t_grid[j]) - u[j,1]*exp(-rates[j]*t_grid[j])), # [1, 1]
      matrix(spot*exp(-dividends[j]*t_grid[j]) - u_j[1]*exp(-rates[j]*t_grid[j])), # [1, 1]
      matrix(0) # [1,1]
    )
    
    if(nrow(LHS_t)!=nrow(rhs)){stop("non conformable size")}
  
    x <- try(solve.QP(Dmat = B, dvec = y, Amat = t(LHS_t), bvec = rhs, meq=(nk-2))$solution,
             silent = TRUE)
    
    if(inherits(x, "try-error")) {
      
      lambda = lambda*1e3
      B <- rbind( cbind(diag(weights), matrix(0, nk, nk-2)),
                  cbind(matrix(0, nk-2, nk), lambda*R) )
      warning(paste("lambda set to ", lambda, "at time ", j))
      x <- solve.QP(Dmat = B, dvec = y, Amat = t(LHS_t), bvec = rhs, meq=(nk-2))$solution
      
      #if(inherits(x, "try-error")) {
      #  warning(paste("quadprog failed at time ", j))
      #  x <- rep(NA, 2*nk-2)
      #}
    } 

    # check: is cubic spline?
    # if(!all.equal(t(Q)%*%x[1:nk], R%*%x[(nk+1):(2*nk-2)])) stop('cubic spline check not passed')
    
    # unconst <- solve.QP(Dmat = B, dvec = y, Amat = (A), bvec = matrix(0, nk-2, 1), meq=(nk-2))$solution
    # x <- unconst
    # 100*(p_j/x[1:nk]-1)
    # plot(diffMat(u_j,2)$mat%*%p_j,type='l')
    # lines(c(0,x[(nk+1):(2*nk-2)],0), col='blue')
    
    # g[j,] <- x[1:nk]
    # gamma[j,2:(nk-1)] <- x[(nk+1):(2*nk-2)]
    g[j,sig.space[j,]] <- x[1:nk]
    gamma[j,sig.space[j,]] <- c(0, x[(nk+1):(2*nk-2)], 0)
    
    # Set small values to zero 
    # gamma[j,which(gamma[j,1:nk]<2.2e-16)] <- 0
    gamma[j,which(gamma[j,]<2.2e-16)] <- 0
    g[j,which(g[j,]<2.2e-16)] <- 0
    
    # CONDITIONS:  
    # ((cbind(matrix(0, nk-2, nk), diag(rep(1,nk-2))))%*%x>=matrix(0, nk-2, 1)) # 1
    # (c(-1/h[1], 1/h[1], rep(0,nk-2), -h[1]/6, rep(0,nk-3))%*%x)>=matrix(-exp(-rates[j]*t_grid[j])) # 2
    # ((c(rep(0,nk-2), 1/h[nk-1], -1/h[nk-1], rep(0,nk-3), -h[nk-1]/6))%*%x)>=0 # 3 
    # nocs_lhs%*%x>=matrix(nocs_rhs) # 4
    # (c(1, rep(0, 2*nk-3))%*%x)>=matrix(spot*exp(-dividends[j]*t_grid[j]) - u_j[1]*exp(-rates[j]*t_grid[j])) # 5 - FALSE
    # c(rep(0, nk-1), 1, rep(0, nk-2))%*%x>=0 # 6 
  }
  
  return(list(p = g, gamma = gamma, K = u))

}

