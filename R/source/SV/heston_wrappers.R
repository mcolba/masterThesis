
###  heston_wrappers.R  ###

dyn.load('D:/Dropbox/Thesis/Code/hestonCalibrator_bis/x64/Release/hestonCalibrator.dll')
# dyn.load('D:/Dropbox/Thesis/Code/hestonCalibrator_bis/x64/Debug/hestonCalibrator.dll')


hestonPricer <- function(par, S, K, mat, r, q) {  
  #
  #   par = [5]-vector of heston parametes 
  #   spot 
  #   k = [n]-vector of strikes 
  #   mat = [n]-vector of maturities 
  #   r = [n]-vector of continously conpounded rates 
  #   q = [n]-vector of continously conpounded dividend yields 
  #
  
  n <- length(K)
  m <- length(par)
  
  # if yields is a constant make it a vector  
  if(length(r)==1){r=rep(r,n)}
  if(length(q)==1){r=rep(q,n)}
  
  # check dim 
  if(n!=length(r)|n!=length(q)|n!=length(K)|n!=length(mat)) stop('wrong dim')
  if(length(par)!=5) warning('wrong parameter dim')
  
  # compute adjusted spot  
  spot_adj <- S*exp(-q*mat)
  
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
  #   par = [5]-vector of heston parametes 
  #   spot 
  #   k = [n]-vector of strikes 
  #   mat = [n]-vector of maturities 
  #   r = [n]-vector of continously conpounded rates 
  #   q = [n]-vector of continously conpounded dividend yields 
  #
  
  n <- length(K)
  m <- length(par)
  
  # if yields is a constant make it a vector  
  if(length(r)==1){r=rep(r,n)}
  if(length(q)==1){r=rep(q,n)}
  
  # check dim 
  if(n!=length(r)|n!=length(q)|n!=length(K)|n!=length(mat)) stop('wrong dim')
  if(length(par)!=5) warning('wrong parameter dim')
  
  # compute adjusted spot  
  spot_adj <- S*exp(-q*mat)
  
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
  #   par = [5]-vector of heston parametes 
  #   spot 
  #   k = [n]-vector of strikes 
  #   mat = [n]-vector of maturities 
  #   r = [n]-vector of continously conpounded rates 
  #   q = [n]-vector of continously conpounded dividend yields 
  #
  
  n <- length(K)
  m <- length(par)
  
  # if yields is a constant make it a vector  
  if(length(r)==1){r=rep(r,n)}
  if(length(q)==1){r=rep(q,n)}
  
  # check dim 
  if(n!=length(r)|n!=length(q)|n!=length(K)|n!=length(mat)) stop('wrong dim')
  if(length(par)!=5) warning('wrong parameter dim')
  
  # compute adjusted spot  
  spot_adj <- S*exp(-q*mat)
  
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
  #
  #  stopCriteria: 
  #   6 -> 'Solved: stopped by small ||e||_2'
  #   1 -> 'Solved: stopped by small gradient t(J)e'
  #   2 -> 'Solved: stopped by small change Dp'
  
  n <- length(K)
  m <- length(guess)
  
  # if yields is a constant make it a vector  
  if(length(r)==1){r=rep(r,n)}
  if(length(q)==1){r=rep(q,n)}
  
  # check dim 
  if(n!=length(price)|n!=length(r)|n!=length(q)|n!=length(K)|n!=length(mat)) stop('wrong dim')
  if(length(guess)!=5) warning('wrong parameter dim')
  
  # compute adjusted spot  
  spot_adj <- S*exp(-q*mat)
  # create info vec 
  statistics <- rep(0, 7)
  
  # Call c++ function 
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

  # stopping criteria 
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
  # parameters vec 
  par <- fit$par
  names(par) = c('k', 'v_inf', 'sigma', 'rho', 'v_0')
  # errors info vec 
  residuals.info <- fit$statistics[1:4]
  names(residuals.info) = c('||e0||_2','||e*||_2','||t(J)e||_inf','||Dp||_2')
  # speed and # loops info vec 
  time.info <- fit$statistics[5:6]
  names(time.info) = c('CPU.time','iterations')
  
  return(list(par=par,
              residuals.info=residuals.info,
              time.info=time.info,
              stopCriteria=stopCriteria))
}

# dyn.unload('D:/Dropbox/Thesis/Code/hestonCalibrator_bis/x64/Release/hestonCalibrator.dll')
# dyn.unload('D:/Dropbox/Thesis/Code/hestonCalibrator_bis/x64/Debug/hestonCalibrator.dll')

test_heston <- function(){
  
  S = 100
  K = 130 
  mat = 3
  r = 0.5
  q = 0.2
  # Heston parameters 
  k = 0.2
  v_inf = 0.5
  sigma = 0.25
  rho = -0.75 
  v0 = 0.2
  par <- c(k,v_inf,sigma,rho,v0)
  # price 
  c <- hestonPricer(par, S, K, mat, r, q)
  # closed form delta 
  delta.an <- hestonDelta(par, S, K, mat, r, q)
  # numerical delta 
  eps <- 1e-6
  diff <- hestonPricer(par, S+eps, K, mat, r, q)-hestonPricer(par, S-eps, K, mat, r, q)
  delta.num <- diff/(2*eps)
  # compare 
  all.equal(delta.an,delta.num, tolerance=1e-4)
}
