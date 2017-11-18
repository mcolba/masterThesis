
###  Dupire94.R  ### 

source("./source/numericalMethods/fdm.R") 

dupireLV_p <- function(c, k, r, q, c_t, c_k, c_kk, t=NULL){
  #
  # rowas time, columns strikes 
  #
  
  if(!is.vector(k)) stop("wrong dimension")
  
  # if missing, compute derivatives numerically 
  if(missing(c_t) & !is.null(t))
    c_t <- diffMat(t,1)$mat%*%c
  if(missing(c_k))
    c_k <- c%*%t(diffMat(k,1)$mat) 
  if(missing(c_kk))
    c_kk <- c%*%t(diffMat(k,2)$mat)
  
  # dupire formula 
  nom <- c_t + q*c + (r-q)*c_k%*%diag(k)
  denom <- c_kk%*%diag(k*k)
  denom[denom==0] <- NA
  return(sqrt(2*nom/denom))
  
}

dupireLV_iv <- function(iv, k, r, q, iv_t, iv_k, iv_kk, t, s){
  
  # check dimansion 
  if(nrow(iv)!=length(t) | ncol(iv)!=length(k)) stop("wrong dimension")
  
  # if missing, compute derivatives numerically 
  if(missing(iv_t))
    iv_t <- diffMat(t,1)$mat%*%iv
  if(missing(iv_k))
    iv_k <- iv%*%t(diffMat(k,1)$mat)
  if(missing(iv_kk))
    iv_kk <- iv%*%t(diffMat(k,2)$mat)
  
  moneyness = log(s/k)
  std <- sqrt(t) * iv
  v1 <- (r - q + 0.5*iv^2)*t
  d1 <- t(moneyness + t(v1))/std
  kt <- matrix(sqrt(t))%*%k
  nom <- 2*iv_t + iv/t + 2*(r-q)*iv_k%*%diag(k)
  denom <- (iv_kk - d1*sqrt(t)*iv_k^2 + ((1/kt + d1*iv_k)^2)/iv)%*%diag(k*k)
  denom[denom==0] <- NA
  
  return(sqrt(nom/denom))

}

test_dupireLV <- function(){
  
  source("./source/BSformulas.R")
  
  # Data from Andersen et al (1997)
  S0 <- 590
  r <- rep(0.06,10)
  q <- rep(0.0262,10)
  
  # Implied volatility: Rows are maturities (same as Andersen97 p. 23)
  IV = c(
    190,	168,	133,	113,	102,	97,	  120,	142,	169,	200,
    177,	155,	138,	125,	109,	103,  100,  114,	130,	150,
    172,	157,	144,	133,	118,	104,	100,	101,	108,	124,
    171,	159,	149,	137,	127,	113,	106,	103,	100,	110,
    171,	159,	150,	138,	128,	115,	107,	103,	99,	  108,
    169,	160,	151,	142,	133,	124,	119,	113,	107,	102,
    169,	161,	153,	145,	137,	130,	126,	119,	115,	111,
    168,	161,	155,	149,	143,	137,	133,	128,	124,	123,
    168,	162,	157,	152,	148,	143,	139,	135,	130,	128,
    168,	164,	159,	154,	151,	148,	144,	140,	139,	132)/1000
  IV <- matrix(IV,10,10, byrow = T)
  K = c(501.50, 531.00, 560.50, 590.00, 619.50, 649.00, 678.50, 708.00, 767.00, 826.00)
  mat = c(0.175, 0.425, 0.695, 0.94, 1, 1.5, 2, 3, 4, 5)
  
  # compute prices 
  p <- BScall_price2D(S0, K, r, q, IV, mat) # rows are maturities 
  
  source("./source/rates.R") 
  # spot to forward rates 
  r_f <- spotTOforward(r, mat, 'contTOcont')$f
  q_f <- spotTOforward(q, mat, 'contTOcont')$f
  
  # Dupire price-formula  
  arg1 <- list(c=p, # prices  
               k=K, # strikes 
               r=r_f, # rates 
               q=q_f, # dividends 
               t=mat) # expiration 
  LV_p <- do.call(dupireLV_p, arg1)
  
  # Dupire IV-formula
  arg2 <- list(iv=IV, # prices  
               k=K, 
               r=r_f, 
               q=q_f, 
               t=mat, 
               s=S0) # expiration 
  LV_iv <- do.call(dupireLV_iv, arg2)
  
  100*LV_iv/LV_p-100
  (LV_p - LV_iv)*100
  
  source("./source/plotter.R") 
  
  visualizeSurface(mat, K, LV_p)
  open3d()
  visualizeSurface(mat, K, LV_iv)
  
}


