
###  test_AndersenBrothertonRatcliffe.R  ### 

# UNIT TESTING: replicating Andersen & Rupert Brotherton-Ratcliffe (1998)

if(!require("testthat")) install.packages("testthat") # for Solve.tridiag 
if(!require("fields")) install.packages("fields") # for Tps
source("./source/BSformulas.R")
source("./source/LV/mesh.R")
source("./source/LV/Andersen_BrothertonRatcliffe.R")
source("./source/numericalMethods/fdm.R")

context("Andersen 1997 fdm")

# parameters 
S0 = 590
r = 0.06
q = 0.0262
# Implied volatility: Rows are maturities (Andersen97 p. 23)
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
# option maturity 
exp <- 2 

# compute prices 
p_initial <- BScall_price2D(S0, K, r, q, IV, mat) # rows are maturities 

# construct grid (equally spaced in log strikes)
t_grid <- seq(mat[1], exp, length.out = 25) 
k_grid <- exp(seq(5.276, 7.553, length.out = 65))
s.idx <- 32
k_grid[s.idx] <- S0 # make S0 a point of the grid, beta = 32 
IVbig <- matrix(NA,length(k_grid),length(t_grid)) # IVS matrix 

# Significant Space 
# Recall, 2SD = 95%, 3SD = 99.7%, 4SD = 99.99%
sl <- S0*exp((r-q-0.5*0.2^2)*t_grid - 3.1*0.2*sqrt(t_grid)) 
su <- S0*exp((r-q-0.5*0.2^2)*t_grid + 3.1*0.2*sqrt(t_grid)) 

# indexing 
significant <- sapply(1:length(t_grid), FUN=function(j){sl[j]<k_grid&k_grid<su[j]}) # 2D
significant[c(1:2,64:nrow(significant)),] <- FALSE  
internal2D <- sapply(1:length(t_grid),FUN=function(j){K[1]<=k_grid&k_grid<=K[length(K)]})
internal1D <-K[1]<=k_grid & k_grid<=K[length(K)]

# INTRAPOLATION (TPS) - Total variance
xy <- expand.grid(mat, K)
tpsFit <- Tps(xy, as.vector(mat*(IV^2)), lambda=1e-13)
new.xy <- expand.grid(t_grid, k_grid[internal1D])
l <- lapply(1:length(t_grid), FUN=function(j){cbind(t_grid[j],k_grid[internal1D])})
new.xy <- do.call(rbind, l)
IVbig[internal2D] <- sqrt(predict(tpsFit,x=new.xy)/new.xy[,1])

# EXTRAPOLATION: smooth, gradual flattening of the volatility curve! 
IVbig[,] <- t(flateningEtrapol(x=t(IVbig[internal1D,]), grid=log(k_grid), internal=internal1D, 
                               method='exponential', lambda=4))

# compute call prices 
p <- t(BScall_price2D(S0, k_grid, r, q, t(IVbig), t_grid)) 

# finf a, b, and c. Andersen 98 
r=rep(r,length(t_grid))
q=rep(q,length(t_grid))
abr <- AndersenBrothertonR(t(p), s_grid=k_grid, t_grid, s.idx, r, q, theta=0.5, constrain=NULL, 
                           sig.idx=significant, keepBoundaryC=F, spotRates=T, ns.iv=0.2, max.iv=0.4)

v <- 2*abr$a
sigma <- sqrt(v)
b <- abr$b
A_ini <- abr$A_ini

test_that("Andersen & Rupert Brotherton-Ratcliffe (1997)", {
  
  # TEST 1: 0<=A_ini<=1
  expect_false(any(A_ini<=0),info='Arrow-Debreu prices < 0') 
  expect_false(any(A_ini>1),info='Arrow-Debreu prices > 1')
  
  # TEST 2: sum(A_ini_j) = P_j
  expect_equal(apply(A_ini,2,sum),c(1, exp(-r*t_grid)),info='constrain 31a not satisfied ')
  
  # TEST 3: sum(A_ini_j*S) = S_ini*1/(1+q*t) 
  expect_equal(apply(A_ini*k_grid,2,sum), c(k_grid[s.idx], k_grid[s.idx]*exp(-q*t_grid)),info='constrain 31b not satisfied')
  
  # TEST 4: b -> r-q-0.5v
  r_f <- spotTOforward(r, t_grid, 'contTOcont')$f
  q_f <- spotTOforward(q, t_grid, 'contTOcont')$f
  expect_equal(b, (r_f-q_f-0.5*v), tolerance=1e-2, info='b different from its continous time version') # 1% tolerance

  # TEST 5: constrained (-inf,inf) = unconstrained
  abr.cs <- AndersenBrothertonR(t(p), s_grid=k_grid, t_grid, s.idx, r, q, theta=0.5, constrain=c(0.04, 0.4), 
                                sig.idx=significant, keepBoundaryC=F, spotRates=T, ns.iv=0.2, max.iv=0.4)
  sigma.cs <- sqrt(2*abr.cs$a)
  # are constrain satisfied? 
  expect_false(any(sigma.cs<(0.04-1e-8)))
  expect_false(any(sigma.cs>(0.4+1e-8)))
  # same unconstrained vs? 
  expect_equal(abr.cs$v.unconstrained, abr$v.unconstrained, tolerance=1e-7, 
               info='b different from its continous time version')
  
  # TEST 6: relative price error < 0.5% 
  # FDM pricer 
  fdmp <- rep(NA,10)
  for(i in 1:10){
    fdmp[i] <- fdmABR(K[i], s_grid=k_grid, s.idx, t_grid, abc=abr, r, q, theta = 0.5, n.imp = 0, spotRates = TRUE)$p
  }
  relativePriceError <- fdmp/p_initial[7,]-1
  zeros <- rep(0,10)
  expect_equal(relativePriceError, zeros, tolerance=5e-3)
})

