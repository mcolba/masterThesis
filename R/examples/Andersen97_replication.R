
###  Andersen97_replication.R  ###

# Replication of the paper 'The equity option volatility smile: an implicit finite-difference approach'
# tabe 2 -> price 2 years call option 

rm(list = ls())

source("./source/plotter.R") 
source("./source/numerical/fdm.R")
source("./source/BSformulas.R")
source("./source/LV/Andersen97.R")
source("./source/LV/mesh.R")
source("./source/LV/extrapolation.R")

S0 = 590
r = 0.0600
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

# compute prices 
p_initial <- BScall_price2D(S0, K, r, q, IV, mat) # rows are maturities 

# option maturity 
exp <- 2 

# construct grid
t_grid <- seq(mat[1], exp, length.out = 25) 
k_grid <- exp(seq(5.276, 7.553, length.out = 65))
s.idx <- 32

k_grid[s.idx] <- S0 # make S0 a point of the grid, beta = 32 
# IVS matrix 
IVbig <- matrix(NA,length(k_grid),length(t_grid))

# Significant Space 
# Recall 2SD = 95%, 3SD = 99.7%, 4SD = 99.99%
sl <- S0*exp((r-q-0.5*0.2^2)*t_grid - 3.1*0.2*sqrt(t_grid)) 
su <- S0*exp((r-q-0.5*0.2^2)*t_grid + 3.1*0.2*sqrt(t_grid)) 

# indexing 
significant <- sapply(1:length(t_grid), FUN=function(j){sl[j]<k_grid&k_grid<su[j]}) # 2D
significant[c(1:2,(nrow(significant)-1):nrow(significant)),] <- FALSE  
internal2D <- sapply(1:length(t_grid),FUN=function(j){K[1]<=k_grid&k_grid<=K[length(K)]})
internal1D <-K[1]<=k_grid&k_grid<=K[length(K)]

# INTRAPOLATION (TPS) - Total variance 
xy <- expand.grid(mat, K)
tpsFit <- Tps(xy, as.vector(mat*(IV^2)), lambda=1e-14, m=3) 
new.xy <- expand.grid(t_grid, k_grid[internal1D])
l <- lapply(1:length(t_grid), FUN=function(j){cbind(t_grid[j],k_grid[internal1D])})
new.xy <- do.call(rbind, l)
IVbig[internal2D] <- sqrt(predict(tpsFit,x=new.xy)/new.xy[,1])

### bicubic -> similar results!!!
# tmpIV <- matrix(NA,length(t_grid),length(K))
# for(i in 1:length(K)){
#   tmpIV[,i] <- spline(mat, IV[,i], method = "natural", xout=t_grid)$y
# }
# IVbig[internal1D,] <- NA 
# for(i in 1:length(t_grid)){
#   IVbig[internal1D,i] <- spline(K, tmpIV[i,], method = "natural", 
#                                 xout=k_grid[internal1D])$y
# }

# EXTRAPOLATION 
# smooth, gradual flattening of the volatility curve
IVbig[,] <- t(flateningEtrapol(x=t(IVbig[internal1D,]), grid=log(k_grid), internal=internal1D, 
                               method='exponential', lambda=5))
# Figure 2 p. 25
open3d()
visualizeSurface(log(k_grid),t_grid,IVbig, title='Figure 2') 
xy <- expand.grid(log(K),mat[1:7]) 
xyz <- cbind(xy, as.vector(t(IV[1:7,])))
points3d(xyz, size=6) 
aspect3d(2,1,1)

# total variance surface 
open3d()
visualizeSurface(t_grid, log(k_grid), t_grid*t(IVbig^2))

# compute call prices 
p <- t(BScall_price2D(S0, k_grid, r, q, t(IVbig), t_grid)) 

# find a, b, and c. Andersen 98 
r=rep(r,length(t_grid))
q=rep(q,length(t_grid))
constrain = NULL
# constrain=c(0.04, 0.5)
abr <- AndersenBrothertonR(p, s_grid=k_grid, t_grid, s.idx, r, q, theta=0.5, constrain, sig.idx=significant, 
                           imposeBoundryValues='greed_boundry', spotRates=T, ns.iv=0.2, max.iv=0.5)
sigma <- sqrt(2*abr$a)
b <- abr$b
A_ini <- abr$A_ini

# figure 3 p. 25 
open3d()
visualizeSurface(log(k_grid), c(0,t_grid), A_ini, title='Figure 3')
# figure 4 p. 27
open3d()
visualizeSurface(log(k_grid),t_grid,sigma, title='Figure 4')
aspect3d(2,1,1)

# FDM pricer 
fdmp <- rep(NA,10)
for(i in 1:10){
  fdmp[i] <- fdmABR(K[i], s_grid=k_grid, s.idx, t_grid, abc=abr, r, q, theta = 0.5, 
                    n.imp = 0, spotRates = TRUE)$p
}

# TABLE 1 
cbind(
  K = K/S0,
  IV = IV[7,],
  BSth = p_initial[7,],
  fdmp = fdmp,
  diff = round((fdmp-p_initial[7,])*100,2),
  relDiff = round((fdmp/p_initial[7,]-1)*100,2)
)

# camparison with Dupire 
source("./source/LV/Dupire94.R")

# partial derivatives 
p_k <- t(diffMat(k_grid,1)$mat%*%(p))
p_kk <- t(diffMat(k_grid,2)$mat%*%(p))
p_t <- (diffMat(t_grid,1)$mat%*%t(p))

r_f <- spotTOforward(r, t_grid, 'contTOcont')$f
q_f <- spotTOforward(q, t_grid, 'contTOcont')$f
nom <- t(p_t + q_f*t(p) + (r_f-q_f)*p_k%*%diag(k_grid))
v.dup <- t(dupireLV_p(t(p), k_grid, r_f, q_f, p_t, p_k, p_kk))

sigma.dup <- sqrt(v.dup)
sigma.dup[sigma.dup>0.5] <- 0.5
sigma.unc <- sqrt(abr$v.unconstrained*2)
sigma.unc[sigma.unc>0.5] <- 0.5

# plot dupire local vol
open3d()
layout3d(matrix(1:2,1,2) , sharedMouse = TRUE)
visualizeSurface(log(k_grid), t_grid, sigma.dup, title='dupire')
aspect3d(2,1,1)
next3d()
visualizeSurface(log(k_grid), t_grid, sigma.unc, title='Andersen')
aspect3d(2,1,1)


# CONSTRAINED vs unconstrained 
abr.cs <- AndersenBrothertonR(p, s_grid=k_grid, t_grid, s.idx, r, q, theta=0.5, 
                              constrain=c(0.04, 0.4), sig.idx=significant, 
                              imposeBoundryValues='greed_boundry', spotRates=T, 
                              ns.iv=0.2, max.iv=0.4)
sigma.cs <- sqrt(2*abr.cs$a)
# plot 
open3d()
layout3d(matrix(1:2,1,2) , sharedMouse = TRUE)
visualizeSurface(log(k_grid),t_grid,sigma.cs, title='constrained')
aspect3d(2,1,1)
next3d()
visualizeSurface(log(k_grid),t_grid,sigma, title='unconstrained')
aspect3d(2,1,1)


