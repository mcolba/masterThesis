
#  Andersen97.R 


# load/install packages 
if(!require("limSolve")) install.packages("limSolve") # for Solve.tridiag 
if(!require("quadprog")) install.packages("quadprog") # for quadratic program. 
# include R files 
source("./source/numerical/fdm.R") 


AndersenBrothertonR <- function(p, s_grid, t_grid, s.idx, r, q, theta, constrain=NULL, 
                                sig.idx=NULL, imposeBoundryValues='greed_boundry', 
                                spotRates=T, ns.iv=NULL, max.iv=NULL){
  # 
  # The function produces the local volatility surface via the implicit finite 
  # difference approach proposed by Andersen et al. (1997). The grid not need 
  # to be equally-spaced. 
  #
  # ARGUMENTS: 
  # * p = [NxM]-matrix of cal option prices 
  # * s_grid = [N]-vector of grid points in the strike-direction. The initial 
  #            price (S0) must be a point of the grid.  
  # * t_grid = [M]-vector of grid points in the time-direction 
  # * s.idx = indexing of the initial price, such that s_grid[s.idx]=S0.  
  # * r = [M]-vector of interest rates  
  # * q = [M]-vector of dividend yields
  # * theta = parameter controlling the finite difference scheme. if 0  -> fully implicit,
  #           if 1 -> fully explicit, if 0.5 Crank-Nicolson. 
  # * constrain = [2]-vector containing min e max IV values to impose in the cinstrained 
  #               method. Default is NULL, which performs the unconstrained method. 
  # * sig.idx = [NxM]-matrix of boolean values indicating the significant grid space.  
  # * imposeBoundryValues = see Andersen et al. (1997)
  # * spotRates = boolean indicating weather r and q are spot or forward rates. 
  # * ns.iv = fixed IV value for the non-significant space.  
  # * max.iv = cap value for the IV.  
  # 
   
  # dim
  N <- length(s_grid)
  M <- length(t_grid)
  
  # if significant space is not provided use all grid excepd 1st and last rows 
  if(missing(sig.idx)){
    sig.idx <- rbind(rep(FALSE,M), matrix(TRUE,N-2,M), rep(FALSE,M))
  } else {
    if(!identical(dim(sig.idx), dim(p))) stop('p and sig.idx have different dim')
  }  
  
  # impose the right ordering: increasing s_grid and t_grid
  idxS <- order(s_grid)  
  s_grid <- s_grid[idxS]  
  idxT <- order(t_grid)
  t_grid <- t_grid[idxT] 
  p <- p[idxS,idxT]; sig.idx <- sig.idx[idxS,idxT] 
  r <- r[idxT]; q <- q[idxT]
  
  # Trasform from/to forward retes 
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
  
  # intervals 
  dt <- c(t_grid[1], t_grid[2:M] - t_grid[1:(M-1)])
  # dx <- x_grid[2:N] - x_grid[1:(N-1)] # dx_i = x_(i+1) - x_i 
  # dx_i <- dx[2:(N-1)]; dx_im1 <- dx[1:(N-2)]
  ds <- s_grid[2:N] - s_grid[1:(N-1)]
  ds_i <- ds[2:(N-1)]; ds_im1 <- ds[1:(N-2)]
  # diff Matrices 
  D1 <- diffMat(x_grid, 1)$mat
  D2 <- diffMat(x_grid, 2)$mat
  
  # Arrow-Debreu securities
  A_ini <- matrix(NA,N,M+1)
  
  # NB: by default call values at the boundries are imposed (see Anderson 97 end note 9)  
  if(imposeBoundryValues=='greed_boundry'){
    
    # inpose values on the upper and lower x boundries  
    p[1,] <- s_grid[s.idx]*exp(-q_s*t_grid)-s_grid[1]*exp(-r_s*t_grid)
    p[N,] <- 0
    
  } else if(imposeBoundryValues=='significant_boundry'){
    
    # impose values for all non significant points 
    for(i in 1:M){
      for(j in 1:N){
        if(!sig.idx[j,i] & j<s.idx){
          p[j,i] <- s_grid[s.idx]*exp(-q_s[i]*t_grid[i])-s_grid[j]*exp(-r_s[i]*t_grid[i])
        }
        if(!sig.idx[j,i] & j>s.idx){
          p[j,i] <- 0
        }
      }
    }
  }

  # Compute A_ini internal As 
  A_ini[-c(1,N),-1] <- (ds_im1*p[3:N,] - (ds_im1+ds_i)*p[2:(N-1),] + ds_i*p[1:(N-2),])/
    (ds_im1*ds_i)
  # Compute A_ini boundry As 
  A_ini[1,-1] <- (s_grid[2]*exp(-r_s*t_grid) - s_grid[s.idx]*exp(-q_s*t_grid) + p[2,])/ds[1]
  A_ini[N,-1] <- p[N-1,]/ds[N-1]   
  A_ini[s.idx,1] <- 1; A_ini[-s.idx,1] <- 0 
  
  # Discrete transition densities
  pj <- c(1,exp(-r_s*t_grid))
  den <- 0.5*c(ds[1], (ds_im1+ds_i), ds[N-1])  
  adPseudoProb <- (A_ini%*%diag(1/pj))

  # chek for very low/negative pseudo probabilities 
  lowADpp <- sum(adPseudoProb[sig.idx]<1e-15)
  if(lowADpp!=0){ 
    warning(paste(lowADpp, 'Grid points with Arrow-Debreu pseudo probabilities < 1e-15 were used'))
  } 
  
  #quantile(exp(-r_s*t_grid) - apply(A_ini,2,sum)[-1])
  #quantile(s_grid[s.idx]*exp(-q_s*t_grid) - apply(A_ini*s_grid,2,sum)[-1])
  
  # c (value term): adjusted to fit all bond prices
  c <- - r_f 
  
  # b (first order derivative term): it's adjusted to fit all forward prices 
  fqDF <- 1/(1+q_f*dt) # dividend forwad discount factor  
  const <- dt*((1-theta)*fqDF + theta) # constant in the x-space  
  # 1st and 2nd derivaives of the s_grid wrt x  
  # s_x <- D1%*%s_grid
  # s_xx <- D2%*%s_grid
  s_x <- s_xx <- s_grid # analitical solution 
  # coefficients 
  a_tilde <- matrix(s_grid/s_x)%*%((fqDF*(1+r_f*dt)-1)/const)
  b_tilde <- -0.5*diag(as.vector(s_xx/s_x))
  
  # a (second order derivative term)  
  v <- matrix(NA,N,M)
  v.unconstrained <- matrix(NA,N,M)
  
  for(j in 1:M){
    
    A_hat <- (1-theta)*A_ini[,j+1] + theta*A_ini[,j]
    Lambda <- dt[j]*t(D1)%*%diag(a_tilde[,j])
    Psi <- dt[j]*(0.5*t(D2) + t(D1)%*%b_tilde)%*%diag(A_hat)
    Tau <- matrix((1+r_f[j]*dt[j])*A_ini[,j+1]-A_ini[,j]) - Lambda%*%matrix(A_hat)
    
    # Reduce dim to the significant space!
    Psi_small <- Psi[sig.idx[,j],sig.idx[,j]]
    Tau_small <- Tau[sig.idx[,j]]
    n <- length(Tau_small)
    
    # LV to ns.iv out of the significant space!!!  
    if(!is.null(ns.iv)){
      low <- which(c(sig.idx[1:(N-1),j]==FALSE&sig.idx[2:N,j]==TRUE, FALSE))
      upp <- which(c(FALSE, (sig.idx[1:(N-1),j]==TRUE&sig.idx[2:N,j]==FALSE)))
      # Tau_small <- Tau_small - ((ns.iv^2)*Psi[sig.idx[,j],low] + (ns.iv^2)*Psi[sig.idx[,j],upp]) 
      Tau_small <- Tau_small - ((0.2^2)*Psi[sig.idx[,j],low] + (0.2^2)*Psi[sig.idx[,j],upp])   #### ?????????????????????    
    }
    
    if(is.null(constrain)){
      
      lhs <- Psi_small
      rhs <- Tau_small
      v[sig.idx[,j],j] <- Solve.tridiag(lhs[cbind(2:(n),1:(n-1))], lhs[cbind(1:(n),1:(n))], 
                                        lhs[cbind(1:(n-1),2:(n))], rhs) 
      v.unconstrained[sig.idx[,j],j] <- v[sig.idx[,j],j]
      
    } else if(is.numeric(constrain) && length(constrain)==2){
      
      v_min <- rep(constrain[1]^2,n); v_max <- rep(constrain[2]^2,n) 
      
      # Constrain on min and max v, solve quadratic program (Phi = I):
      # min. 0.5 y_j' Q y_j + c_j' y_j
      # s.t. 0 <= y_j <= y_max
      # where: y_j <- v_j - v_min
      
      Q <- t(Psi_small)%*%Psi_small
      cvec <- t(Psi_small)%*%(-Tau_small+Psi_small%*%v_min)
      y_max <- v_max-v_min
      
      A <- t(rbind(diag(-1,n),diag(1,n)))
      bvec <- c(-y_max, rep(0,n))
      
      y_j <- solve.QP(Dmat=Q, dvec=-cvec, Amat=A, bvec=bvec)
      # constrained solution 
      v[sig.idx[,j],j] <- (y_j$solution + v_min)
      # non-constrained solution
      v.unconstrained[sig.idx[,j],j] <- (y_j$unconstrained.solution + v_min)
      
    } else { stop('constrain has wrong dimension') }
  }
  
  # cap high variance 
  if(!is.null(max.iv)){
    nExtreme <- sum(!is.na(v) & v>max.iv)
    if(nExtreme!=0){
      v[!is.na(v) & v>max.iv] <- max.iv
      warning(paste(nExtreme, 'variances exceded the cap value of', max.iv))
    }    
  }
  
  # set negative variances to NA 
  nNeg <- sum(!is.na(v) & v<0)
  if(nNeg!=0){
    v[!is.na(v) & v<0] <- NA
    warning(paste(nNeg, 'negative variances were set to NA'))
  }
  
  # set NAs to non significant IV
  v[is.na(v)] <- ns.iv^2
  
  # compute a and b 
  a <- 0.5*v
  b <- a_tilde + b_tilde%*%v
  
  return(list(a=a, b=b, c=c, A_ini=A_ini, adPseudoProb=adPseudoProb, v.unconstrained=v.unconstrained))
  
}

