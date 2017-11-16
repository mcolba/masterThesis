
###  fdm.R  ###

if (!require("SoDA")) install.packages("SoDA") 
if (!require("limSolve")) install.packages("limSolve") # for Solve.tridiag 
source("./source/rates.R") # for spotTOforward

fdmABR <- function(K, s_grid, s.idx, t_grid, abc, r, q, theta = 0.5, n.imp = 0, 
                   spotRates = TRUE){
  #
  # The function perform a finite difference method to price a call option trough the
  # backward PDE. 
  #
  # ARGUMENTS:
  # * s.idx = index of the s_grid corresponding to the desired spot,                    
  #   i.e. s_grid[s.idx] = spot.  
  # * K = strike 
  # * s_grid = [1 x N]-vector, the spot-dimension mesh.  
  # * t_grid = [1 x M]-vector, the t-dimension mesh.
  # * p = [M x N]-matrix of call option prices.
  # * r,q = [1 x M]-vectos of zero rates and dividend yields. 
  # * theta = 0.5 for the Crank-Nicolson method or 1 for the implicit method, 
  #   default is 0.5 
  # * n.imp = number of initial filly implicit steps, default is 0. 
  # 
  # VALUE: 
  # * the function returns a list containing the following objects: 
  #   * $p = price of the call option with spot = s_grid[s.idx]
  #   * $v = vector of all prices 
  #   * $delta  = vector of all deltas 
  # 
  
  # dim
  N <- length(s_grid)
  M <- length(t_grid)

  # theta vector 
  theta <- c(rep(1,n.imp), rep(theta,M-n.imp))
  
  # impose the right ordering
  idxS <- order(s_grid)
  s_grid <- s_grid[idxS]
  idxT <- order(t_grid)
  t_grid <- t_grid[idxT]
  # p <- p[idxS,idxT]
  r <- r[idxT]
  q <- q[idxT]
  
  # Trasform to forward retes 
  if(spotRates==TRUE){
    r_s <- r; q_s <- q
    r_f <- spotTOforward(r_s, t_grid, type = 'contTOsimp')$f
    q_f <- spotTOforward(q_s, t_grid, type = 'contTOsimp')$f
  } else {
    r_f <- r; q_f <- q
    r_s <- forwardTOspot(r_f, t_grid, type = 'simpTOcont')$s
    q_s <- forwardTOspot(q_f, t_grid, type = 'simpTOcont')$s
  }
  
  # change of variable S -> ln(S)
  x_grid <- log(s_grid)
  
  dt <- c(t_grid[1], t_grid[2:M] - t_grid[1:(M-1)])
  # differenciation matrices 
  D1 <- diffMat(x_grid,1)$mat
  D2 <- diffMat(x_grid,2)$mat
  D1.s <- diffMat(s_grid,1)$mat
  
  # extract pde's diffusion, convection, and reaction coefficients
  a <- abc$a
  b <- abc$b
  c <- abc$c
  
  # boundries (Dirichlet) 
  lowerBound <- rep(0,M+1)
  upperBound <- c(s_grid[N]-K, s_grid[N]*exp(-q_s*t_grid) - K*exp(-r_s*t_grid)) 
  # or??? 
  upperBound <- c((s_grid[N]*exp(-q_s*t_grid)-K*exp(-r_s*t_grid))[M:1], s_grid[N]-K) 
  
  # initial condition (call payoff at maturity)
  v <- matrix(pmax(s_grid[-c(1,N)]-K,0))

  for (j in M:1){
    
    # the system of N ODEs is: 
    # ((1+dt*r)I - (1-theta)*L_j)*H_j = (theta*L_j + I)*H_(j+1) + B_j
    # where L_j = dt((A_j)D2 + (B_j)D1), A = diag(a_j). 
    
    L <- dt[j]*(diag(a[,j])%*%D2 + diag(b[,j])%*%D1)
    B <- L[-c(1,N),1]*((1-theta[j])*lowerBound[j] + theta[j]*lowerBound[j+1]) + 
      L[-c(1,N),N]*((1-theta[j])*upperBound[j] + theta[j]*upperBound[j+1])
    I <- diag(1,(N-2),(N-2))
    
    lhs <- ((1-c[j]*dt[j])*I - (1-theta[j])*L[-c(1,N),-c(1,N)])
    rhs <- (theta[j]*L[-c(1,N),-c(1,N)] + I)%*%v + B
    
    if(!theta[j]==0){
      # solve tridiagonal system 
      v <- Solve.tridiag(lhs[cbind(2:(N-2),1:(N-3))], lhs[cbind(1:(N-2),1:(N-2))], 
                         lhs[cbind(1:(N-3),2:(N-2))], rhs)
    } else {v <- B}
  }
  
  # add boundry values 
  v <- c(lowerBound[1], v, upperBound[1])
  # compute delda 
  delta <- (D1.s%*%v)[s.idx]
  
  # return 
  return(list(p=v[s.idx], v = v, delta = delta, s_grid = s_grid))
  
} 

fdmBack <- function(K, s_grid, s.idx, t_grid, lv, r, q, theta = 0.5, n.imp = 0, 
                    ab.adjust = TRUE, spotRates = FALSE){
  #
  # The function perform a finite difference method to price a call option trough the
  # backward PDE. 
  #
  # ARGUMENTS:
  # * s.idx = index of the s_grid corresponding to the desired spot,                    
  #   i.e. s_grid[s.idx] = spot.  
  # * K = strike 
  # * s_grid = [1 x N]-vector, the spot-dimension mesh.  
  # * t_grid = [1 x M]-vector, the t-dimension mesh.
  # * lv = [M x N]-matrix of local volatilities.
  # * r,q = [1 x M]-vectos of zero rates and dividend yields. 
  # * theta = 0.5 for the Crank-Nicolson method or 0 for the implicit method, 
  #   default is 0.5 
  # * n.imp = number of initial filly implicit steps, default is 0. 
  # 
  # VALUE: 
  # * the function returns a list containing the following objects: 
  #   * $p = price of the call option with spot = s_grid[s.idx]
  #   * $v = vector of all prices 
  #   * $delta  = vector of all deltas 
  # 

  # dim
  N <- length(s_grid)
  M <- length(t_grid)
  
  # impose the right ordering
  idxS <- order(s_grid)
  s_grid <- s_grid[idxS]
  idxT <- order(t_grid)
  t_grid <- t_grid[idxT]
  r <- r[idxT]; q <- q[idxT]
  lv <- lv[idxS,idxT]
  
  # theta vector 
  theta <- c(rep(theta,M-n.imp), rep(0,n.imp))

  # Trasform to forward retes 
  if(spotRates==TRUE){
    r_s <- r; q_s <- q
    r_f <- spotTOforward(r_s, t_grid, type = 'contTOsimp')$f
    q_f <- spotTOforward(q_s, t_grid, type = 'contTOsimp')$f
  } else {
    r_f <- r; q_f <- q
    r_s <- forwardTOspot(r_f, t_grid, type = 'simpTOcont')$s
    q_s <- forwardTOspot(q_f, t_grid, type = 'simpTOcont')$s
  }
  
  dt <- c(t_grid[1], t_grid[2:M] - t_grid[1:(M-1)])
  # dx <- s_grid[2:N] - s_grid[1:(N-1)] # dx_i = x_(i+1) - x_i 
  # dx_i <- dx[2:(N-1)]
  # dx_im1 <- dx[1:(N-2)]
  
  # differenciation matrices 
  D1 <- diffMat(s_grid, 1)$mat
  D2 <- diffMat(s_grid, 2)$mat

  # NB: we use the change of variable tao = T - t, i.e. t_grid = 0, ..., T
  # NO!
  
  if(ab.adjust==TRUE){
    
    # c (value term): it's adjusted to fit all bond prices
    c <- -r_f/(1+(1-theta)*r_f*dt)
    
    # b (first order derivative term): it's adjusted to fit all forward prices 
    # NB: D1%*%s_grid = 1, D2%*%s_grid = 0
    fqDF <- 1/(1+q_f*dt) # dividend forwad discount factor  
    v1 <- (fqDF-1)/(dt*(theta + (1-theta)*fqDF))
    b <- matrix(s_grid)%*%t(matrix(v1-c)) 
    # check:  
    # quantile((matrix(s_grid) %*% t(matrix(r_f - q_f))) - b)
    
  } else {
    
    c <- -r_f
    b <- matrix(s_grid)%*%t(matrix(r_f - q_f))
  }

  # a (second order derivative term)  
  a <- 0.5*(lv*lv)*(s_grid*s_grid) 
  
  # boundries (Dirichlet) 
  lowerBound <- rep(0,M+1)
  upperBound <- c(s_grid[N]*exp(-q_s[M:1]*t_grid[M:1]) - K*exp(-r_s[M:1]*t_grid[M:1]), s_grid[N]-K) 
  
  # initial condition (call payoff at maturity)
  v <- matrix(pmax(s_grid[-c(1,N)]-K,0))
  
  # loop: solve backward (from M to 1) 
  for (j in M:1){
    
    # the system of N ODEs is: 
    # (I + theta*dt*L_(j+1))*v_(j+1) = (I - (1-theta)*dt*L_j)*v_j - R_j
    # or, A %*% v_(j+1) = B 
    
    # we carry out the discretisation at t_(j+theta)
    L <- L_next <- diag(a[,j])%*%D2 + diag(b[,j])%*%D1 + c[j]*diag(1,N,N)
    
    # Dirichlet condition 
    R <- (1-theta[j])*dt[j]*(L[-c(1,N),1]*lowerBound[j] + L[-c(1,N),N]*upperBound[j]) + 
      theta[j]*dt[j]*(L_next[-c(1,N),1]*lowerBound[j+1] + L_next[-c(1,N),N]*upperBound[j+1])
    
    LHS <- (diag(1,(N-2),(N-2)) - (1-theta[j])*dt[j]*L[-c(1,N),-c(1,N)])
    RHS <- (diag(1,(N-2),(N-2)) + theta[j]*dt[j]*L_next[-c(1,N),-c(1,N)])%*%v + R
    
    if(!theta[j]==1){
      
      v <- Solve.tridiag(LHS[cbind(2:(N-2),1:(N-3))], LHS[cbind(1:(N-2),1:(N-2))], 
                         LHS[cbind(1:(N-3),2:(N-2))], RHS)
    } else {v <- RHS}
    
    # L <- L_next
  }
  
  # add boundry values 
  v <- c(lowerBound[1], v, upperBound[1])
  
  return(list(p=v[s.idx], v = v, delta = (D1%*%v)[s.idx], s_grid = s_grid))
  
}

diffMat <- function(grid, order = 1, accurate.boundry=FALSE){
  #
  #   given a function y(x), and the respective vectors 
  #   Y = [y_1, y_2, ... , y_n] and X = [x_1, x_2, ... , x_n], 
  #   We have:     
  #   Y %*% diffMat(X,1) = [y'_1, y'_2, ... , y'_n]'   
  #
  # Description:
  #   We use a 3-point central difference everiwere exept for the boundaries,
  #   - 1st and last first order derivatives: one-sided 2-point approx. 
  #   - 1st and last second order derivative: one-sided 3-point approx (the finite
  #     difference at 0 and N+1 is the same as  the central difference at 1 and N).
  #
  # ref: http://www.m-hikari.com/ijma/ijma-password-2009/ijma-password17-20-2009/bhadauriaIJMA17-20-2009.pdf
  
  N <- length(grid)
  if(N<3){ stop('insufficient grid points') }
  dx <- grid[2:N] - grid[1:(N-1)] # dx_i = x_(i+1) - x_i 
  # create diagonals for the INTERIOR POINTS of the differenciation matrices 
  dx_i <- dx[2:(N-1)]
  dx_im1 <- dx[1:(N-2)]
  
  # first order 
  if(order==1){
    # lower 
    d1_l <- -dx_i/(dx_im1*(dx_im1+dx_i))
    # central 
    d1_c <- (dx_i-dx_im1)/(dx_i*dx_im1) 
    # upper 
    d1_u <- dx_im1/(dx_i*(dx_im1+dx_i))
    # create differentiation matrices:
    # add 2p. forward/backward differences to D1's BOUNDRY POINTS 
    D1 <- triDiag(c(-1/dx[1], d1_c, 1/dx[N-1]), 
                  c(1/dx[1], d1_u), 
                  c(d1_l, -1/dx[N-1]))
    
    if(accurate.boundry==TRUE){
      # 4p. forward differences to D1's BOUNDRY POINTS
      H_1 <- sum(dx[1:3]) 
      coef1 <- -((2*dx[1]+dx[2])*H_1 + dx[1]*(dx[1]+dx[2]))/
               (dx[1]*(dx[1]+dx[2])*H_1)
      coef2 <- (dx[1]+dx[2])*H_1/
               (dx[1]*dx[2]*(dx[1]+dx[3]))
      coef3 <- -dx[1]*H_1/
               ((dx[1]+dx[2])*dx[2]*dx[3])
      coef4 <- dx[1]*(dx[1]+dx[2])/
               (H_1*(dx[2]+dx[3])*dx[3])
      D1[1,1:4] <- c(coef1,coef2,coef3,coef4) 
      # 4p. backword differences to D2's BOUNDRY POINTS
      H_n <- sum(dx[(N-3):(N-1)])
      coef1 <- -(dx[N-2]+dx[N-1])*dx[N-1]/
               (dx[N-3]*(dx[N-3]+dx[N-2])*H_n)
      coef2 <- H_n*dx[N-1]/
               (dx[N-3]*dx[N-2]*(dx[N-3]+dx[N-1]))
      coef3 <- -H_n*(dx[N-2]+dx[N-1])/
               ((dx[N-3]+dx[N-2])*dx[N-2]*dx[N-1])
      coef4 <- ((3*H_n-dx[N-3])*dx[N-1]+H_n*dx[N-2])/
               (H_n*(dx[N-2]+dx[N-1])*dx[N-1])
      D1[N,(N-3):N] <- c(coef1,coef2,coef3,coef4) 
    }
    
    return(list(mat=D1, l=d1_l, c=d1_c, u=d1_u)) 
  }
  
  # second order 
  else if(order==2){
    # lower 
    d2_l <- 2/(dx_im1*(dx_im1+dx_i)) 
    # central 
    d2_c <- -2/(dx_i*dx_im1) 
    # upper 
    d2_u <- 2/(dx_i*(dx_im1+dx_i)) 
    # create differentiation matrices:
    D2 <- triDiag(c(NA ,d2_c, NA), c(NA, d2_u), c(d2_l, NA))
    
    if(accurate.boundry==FALSE){
      # 3p. forward/backward differences to D2's BOUNDRY POINTS
      D2[c(1,N),c(1:3,(N-2):N)] <- D2[c(2,(N-1)),c(1:3,(N-2):N)]
    } else {
      # 4p. forward differences to D2's BOUNDRY POINTS
      H_1 <- sum(dx[1:3]) 
      coef1 <- 2*(3*dx[1]+2*dx[2]+dx[3])/
               (dx[1]*(dx[1]+dx[2])*H_1)
      coef2 <- -2*(2*dx[1]+2*dx[2]+dx[3])/
               (dx[1]*dx[2]*(dx[1]+dx[3]))
      coef3 <- 2*(dx[1]+H_1)/
               ((dx[1]+dx[2])*dx[2]*dx[3])
      coef4 <- -2*(2*dx[1]+dx[2])/
               (H_1*(dx[2]+dx[3])*dx[3])
      D2[1,1:4] <- c(coef1,coef2,coef3,coef4) 
      # 4p. backword differences to D2's BOUNDRY POINTS
      H_n <- sum(dx[(N-3):(N-1)])
      coef1 <- -2*(dx[N-2]+2*dx[N-1])/
               (dx[N-3]*(dx[N-3]+dx[N-2])*H_n)
      coef2 <- 2*(H_n+dx[N-1])/
               (dx[N-3]*dx[N-2]*(dx[N-3]+dx[N-1]))
      coef3 <- -2*(2*H_n-dx[N-3])/
               ((dx[N-3]+dx[N-2])*dx[N-2]*dx[N-1])
      coef4 <- 2*(dx[N-3]+2*dx[N-2]+3*dx[N-1])/
               (H_n*(dx[N-2]+dx[N-1])*dx[N-1])
      D2[N,(N-3):N] <- c(coef1,coef2,coef3,coef4) 
     }

    return(list(mat=D2, l=d2_l, c=d2_c, u=d2_u)) 
  }
  
  else stop('order must be 1 or 2')

}




