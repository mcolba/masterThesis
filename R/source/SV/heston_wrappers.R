
###  heston_wrappers.R  ###

dyn.load('D:/Dropbox/Thesis/Code/hestonCalibrator/x64/Release/hestonCalibrator.dll')
# dyn.load('D:/Dropbox/Thesis/Code/hestonCalibrator/x64/Debug/hestonCalibrator.dll')

hestonPricer <- function(par, S, K, mat, r, q) {  
  # 
  #   The routine call the c++ code by Cui et al. (2017) for computing the price 
  #   of an European call option under the Heston model.  
  #
  # ARGUMENTS: 
  #   * par = [5]-vector of heston parametes in the order: k, v_inf, volofvol
  #           rho, and v_0.     
  #   * S = underlying spot price 
  #   * K = [n]-vector of strikes 
  #   * mat = [n]-vector of maturities 
  #   * r = [n]-vector of continously conpounded (spot) rates 
  #   * q = [n]-vector of continously conpounded (spot) dividend yields 
  #
  # VALUE: 
  #   * The function returns a [n]-vector of call option prices. 
  # 
  
  n <- length(K)
  m <- length(par)
  
  # if the yield is a constant make it a vector  
  if(length(r)==1){r=rep(r,n)}
  if(length(q)==1){r=rep(q,n)}
  
  # check dim 
  if(n!=length(r)|n!=length(q)|n!=length(K)|n!=length(mat)) stop('wrong dim')
  if(length(par)!=5) stop('wrong parameter dim')
  
  # compute adjusted spot  
  spot_adj <- S*exp(-q*mat)
  
  # call c++ routine 
  prices <- .C('hestonPricer',
               S_adj = as.double(spot_adj), 
               r = as.double(r), 
               K = as.double(K), 
               mat = as.double(mat), 
               par = as.double(par),
               n = as.integer(n), 
               m = as.integer(m),
               p = as.double(rep(0.0,n)))$p
  
  return(prices)
  
}

hestonDelta <- function(par, S, K, mat, r, q) {  
  # 
  #   The routine computes the delta of an European call option under the Heston model.  
  # 
  # ARGUMENTS: 
  #   * par = [5]-vector of heston parametes in the order: k, v_inf, volofvol
  #           rho, and v_0.     
  #   * S = underlying spot price 
  #   * K = [n]-vector of strikes 
  #   * mat = [n]-vector of maturities 
  #   * r = [n]-vector of continously conpounded (spot) rates 
  #   * q = [n]-vector of continously conpounded (spot) dividend yields 
  #
  # VALUE: 
  #   * The function returns a [n]-vector of call option deltas 
  #
  
  n <- length(K)
  m <- length(par)
  
  # if the yield is a constant make it a vector  
  if(length(r)==1){r=rep(r,n)}
  if(length(q)==1){r=rep(q,n)}
  
  # check dim 
  if(n!=length(r)|n!=length(q)|n!=length(K)|n!=length(mat)) stop('wrong dim')
  if(length(par)!=5) stop('wrong parameter dim')
  
  # compute adjusted spot  
  spot_adj <- S*exp(-q*mat)
  
  # call c++ routine 
  delta <- .C('hestonDelta',
              S_adj = as.double(spot_adj), 
              r = as.double(r),
              q = as.double(q),
              K = as.double(K), 
              mat = as.double(mat), 
              par = as.double(par),
              n = as.integer(n), 
              m = as.integer(m),
              p = as.double(rep(0.0,n)))$p
  
  return(delta)
  
}

hestonJacobian <- function(par, S, K, mat, r, q) {  
  # 
  #   The routine call the c++ code by Cui et al. (2017) for computing the jacobian 
  #   of an European call option under the Heston model.  
  #
  # ARGUMENTS: 
  #   * par = [5]-vector of heston parametes in the order: k, v_inf, volofvol
  #           rho, and v_0.     
  #   * S = underlying spot price 
  #   * K = [n]-vector of strikes 
  #   * mat = [n]-vector of maturities 
  #   * r = [n]-vector of continously conpounded (spot) rates 
  #   * q = [n]-vector of continously conpounded (spot) dividend yields 
  #
  # VALUE: 
  #   * The function returns a [n x 5]-matrix of first order derivatives  
  #
  
  n <- length(K)
  m <- length(par)
  
  # if the yield is a constant make it a vector  
  if(length(r)==1){r=rep(r,n)}
  if(length(q)==1){r=rep(q,n)}
  
  # check dim 
  if(n!=length(r)|n!=length(q)|n!=length(K)|n!=length(mat)) stop('wrong dim')
  if(length(par)!=5) stop('wrong parameter dim')
  
  # compute adjusted spot  
  spot_adj <- S*exp(-q*mat)
  
  # call c++ routine 
  jac <- .C('hestonJac',
               S_adj = as.double(spot_adj), 
               r = as.double(r), 
               K = as.double(K), 
               mat = as.double(mat), 
               par = as.double(par),
               n = as.integer(n), 
               m = as.integer(m),
               jac = as.double(rep(0.0,n*m)))$jac
  
  return(jac)
  
}

hestonCalibrator <- function(price, S, K, mat, r, q, guess, printSummary=TRUE) {  
  
  # 
  #   The routine call the c++ code by Cui et al. (2017) for calibrating the Heston 
  #   model to a set of option prices. 
  #
  # ARGUMENTS: 
  #   * price = [n]-vector of European call option prices  
  #   * S = underlying spot price 
  #   * K = [n]-vector of strikes 
  #   * mat = [n]-vector of maturities 
  #   * r = [n]-vector of continously conpounded (spot) rates 
  #   * q = [n]-vector of continously conpounded (spot) dividend yields 
  #   * guess = [5]-vector containing the initial guess for the heston parametes, 
  #             in the order: k, v_inf, volofvol, rho, and v_0. 
  #   * printSummary = boolean stating whether to print or not the summary info. 
  # VALUE: 
  #   * a list with the following components 
  #     * $par = [5]-vector of heston parametes in the order: k, v_inf, volofvol
  #              rho, and v_0. 
  #     * $residuals.info = [4]-vector containing: ||e0||_2, ||e*||_2, 
  #                         ||t(J)e||_inf, ||Dp||_2
  #     * $time.info  = [2]-vector containing the CPU time and the nr. of iterations. 
  #     * $stopCriteria:  
  #       if == 6 -> Solved: stopped by small ||e||_2
  #       if == 1 -> Solved: stopped by small gradient t(J)e
  #       if == 2 -> Solved: stopped by small change Dp
  
  n <- length(K)
  m <- length(guess)
  
  # if the yield is a constant make it a vector  
  if(length(r)==1){r=rep(r,n)}
  if(length(q)==1){r=rep(q,n)}
  
  # check dim 
  if(n!=length(price)|n!=length(r)|n!=length(q)|n!=length(K)|n!=length(mat)) stop('wrong dim')
  if(length(guess)!=5) warning('wrong parameter dim')
  
  # compute adjusted spot  
  spot_adj <- S*exp(-q*mat)
  # create info vec 
  statistics <- rep(0, 7)
  
  # call c++ routine 
  fit <- .C('hestonCalibrator',
            S_adj = as.double(spot_adj), 
            r = as.double(r), 
            K = as.double(K), 
            mat = as.double(mat), 
            price = as.double(price),
            par = as.double(guess),
            statistics = as.double(statistics),
            n = as.integer(n), 
            m = as.integer(m),
            print = as.logical(printSummary))

  # unsolve error messages  
  stopCriteria <- fit$statistics[7]   
  if (stopCriteria==3){ 
    stop('Unsolved: stopped by itmax')
  } else if(stopCriteria==4) { 
    stop('Unsolved: singular matrix. Restart from current p with increased mu')
  } else if(stopCriteria==5) { 
    stop('Unsolved: no further error reduction is possible. Restart with increased mu')
  } else if(stopCriteria==7) { 
    stop('Unsolved: stopped by invalid values, user error')
  }
  
  # parameters  
  par <- fit$par
  names(par) = c('k', 'v_inf', 'sigma', 'rho', 'v_0')
  
  # errors info  
  residuals.info <- fit$statistics[1:4]
  names(residuals.info) = c('||e0||_2','||e*||_2','||t(J)e||_inf','||Dp||_2')
  
  # CPU time and nr. of loops info  
  time.info <- fit$statistics[5:6]
  names(time.info) = c('CPU.time','iterations')
  
  return(list(par=par,
              residuals.info=residuals.info,
              time.info=time.info,
              stopCriteria=stopCriteria))
}

# dyn.unload('D:/Dropbox/Thesis/Code/hestonCalibrator/x64/Release/hestonCalibrator.dll')
# dyn.unload('D:/Dropbox/Thesis/Code/hestonCalibrator/x64/Debug/hestonCalibrator.dll')

test_heston <- function(){
  
  # mkt parameters 
  S = 100; K = 130; mat = 3; r = 0.5; q = 0.2
  # Heston parameters 
  k = 0.2; v_inf = 0.5; sigma = 0.25; rho = -0.75; v0 = 0.2
  par <- c(k,v_inf,sigma,rho,v0)
  # compute price 
  c <- hestonPricer(par, S, K, mat, r, q)
  # compute closed form delta 
  delta.an <- hestonDelta(par, S, K, mat, r, q)
  # compute numerical delta 
  eps <- 1e-6
  diff <- hestonPricer(par, S+eps, K, mat, r, q)-hestonPricer(par, S-eps, K, mat, r, q)
  delta.num <- diff/(2*eps)
  # compare 
  all.equal(delta.an,delta.num, tolerance=1e-4)
}
