
### LVS_example.R ###    

rm(list = ls())

# include R files 
source("./source/BSformulas.R")
source("./source/LV/Dupire94.R")
source("./source/numerical/rootFinder.R")
source("./source/plotter.R") 
source("./source/LV/cubicSpline.R")
source("./source/LV/Andersen97.R")
source("./source/LV/pre_smoother.R")
source("./source/numerical/fdm.R")
source("./source/LV/extrapolation.R")
source("./source/LV/mesh.R")
source("./source/LV/Fengler09.R")

# load data 
load("./data/data.RData")
N <- length(data) # number of days 
t=250

# unlist varables 
iv <- data[[t]]$iv
S0 <- data[[t]]$spot
K <- data[[t]]$strikes
mat <- data[[t]]$expiries
r.ini <- data[[t]]$rates
q.ini <- data[[t]]$dividends
spotIV <- mean(iv[,6])

# original data dim 
nk.ini <- length(K)
nt.ini <- length(mat)

# reorder
tmp <- order(K) # increasing strikes 
K <- K[tmp]
iv <- iv[,tmp]

# compute mkt prices 
p_mkt <- BScall_price2D(S0,K,r.ini,q.ini,iv,mat)


# 0. GRID -------------------------------------------------------------------------------  
mat.idx = 8                   # <<<<<---- set maturity of the option to be priced 
t.max <- mat[mat.idx]     
nk <- 65                      # <<<<---- set grid points in strike-direction 
nt <- 25  # <<<<---- set grid points in the time-direction 
mesh.sd <- 5                 # <<<<---- set s.d. of the boundry values from S0

# do T-GRID (equally spaced)
t_grid <- seq(mat[1], t.max, length.out = nt) 
# interpolate rates (flat forward)
q <- approx(x = mat, y = q.ini*mat, xout = t_grid)$y/t_grid
r <- approx(x = mat, y = r.ini*mat, xout = t_grid)$y/t_grid

# K-mesh boundry values 
s.min <- S0*exp((r[length(r)]-q[length(q)]-0.5*spotIV^2)*t.max - mesh.sd*spotIV*sqrt(t.max)) 
s.max <- S0*exp((r[length(r)]-q[length(r)]-0.5*spotIV^2)*t.max + mesh.sd*spotIV*sqrt(t.max)) 

# (hyperbolic) k_grid 
tmp <- hyperbolicMesh(s.min, s.max, S0, n=nk, alpha=0.01, x_star=S0)
k_grid <- tmp$x
s.idx <- tmp$idx
# or Equally spaced: 
# k_grid <- seq(s.min, s.max, length.out = nk-1)
# k_grid <- c(k_grid, S0)
# s.idx <- order(k_grid)[nk]
#k_grid <- sort(k_grid)

# Significant Space (recall 3.1SD = 99.9%)
sd_0 = 3.5         # <<<<---- set # of sd.s of the significant area (time 0) 
sd_t.max = 3.5      # <<<<---- set # of sd.s of the significant area (t max)
std <- seq(sd_0,sd_t.max,length.out = length(t_grid))
sl <- S0*exp((r-q-0.5*spotIV^2)*t_grid - std*spotIV*sqrt(t_grid)) 
su <- S0*exp((r-q-0.5*spotIV^2)*t_grid + std*spotIV*sqrt(t_grid)) 
sig <- sapply(1:length(t_grid), FUN=function(j){sl[j]<k_grid&k_grid<su[j]})
sig[c(1:2,(nk-1):nk),] <- FALSE  

# internal space (to interpolate) 
internal2D <- sapply(1:length(t_grid),FUN=function(j){min(K)<=k_grid&k_grid<=max(K)})
internal1D <- min(K)<=k_grid&k_grid<=max(K)

# 1. IV INTERPOLATION/EXTRAPOLATION -----------------------------------------------------
# -- 1.1. Arbitrage Free: Fengler (2007) + Benaim ----------------------------------------

nk.fm = 40 # <<<<---- set number of points in the forward moneyness grid  

# Fengler: TPS presmoothing (default uses 30 equally spaced grid points)  
preSmoothed <- TPSpreSmoother(iv, expiries=mat, strikes=K, spot=S0, rates=r.ini,
                              dividends=q.ini, t_grid.new=t_grid,  alpha = NULL, n=nk.fm) 
smoothP <- preSmoothed$p

# Significant space in the forward moneyness grid 
forward <- S0*exp((r-q)*t_grid)
sl_fm <- sl/forward - 0.1 # lower bound 
su_fm <- su/forward + 0.1 # uooer bound  
fm_grid <- preSmoothed$k_grid
sig_fm <- t(sapply(1:length(t_grid), FUN=function(j){sl_fm[j]<fm_grid&fm_grid<su_fm[j]})) # 2D

# Apply Fengler 07 
arguments <- list(p=smoothP,
                  k_grid = preSmoothed$k_grid,
                  t_grid = t_grid,
                  rates = preSmoothed$dividends,
                  dividends = preSmoothed$dividends,
                  forward = preSmoothed$forward,
                  smooth = 1e-12, # <<<<---- set smoothness parameter 
                  spot = S0, 
                  sig.space=NULL)

arbFreeIV <- do.call(solveQuadprog, arguments)

# if the Fengler routine fails -> run only in the significant space!  
if(any(is.na(arbFreeIV$p))){
  
  arguments <- list(p=smoothP,k_grid = preSmoothed$k_grid,t_grid = t_grid,
                    rates = preSmoothed$dividends,dividends = preSmoothed$dividends,
                    forward = preSmoothed$forward,smooth = 1e-3,spot = S0, 
                    sig.space=sig_fm)
  arbFreeIV <- do.call(solveQuadprog, arguments)
  
  # evaluate at desired strikes 
  p_intra <- delta_intra <- gamma_intra <- matrix(NA,sum(internal1D),length(t_grid))
  for(j in 1:length(t_grid)){
    tmp  <- c.spline_predict(k_grid[internal1D], arbFreeIV$K[j,sig_fm[j,]], arbFreeIV$p[j,sig_fm[j,]], 
                             arbFreeIV$gamma[j,sig_fm[j,]])
    p_intra[,j] <- tmp$g
    delta_intra[,j] <- tmp$delta
    gamma_intra[,j] <- tmp$gamma
  }
  
  p_extra <- matrix(NA,nk,nt)
  
  for(i in 1:nt){
    
    tmp.idx1 <- internal1D & sig[,i]
    tmp.idx2 <- sig[,i][internal1D]
    
    p_extra[,i] <-  benaimEtrapol(t(p_intra)[i,tmp.idx2],k_grid,tmp.idx1, S0, r[i], q[i], t_grid[i], 
                                  c_kk.analytic = gamma_intra[tmp.idx2,i], 
                                  c_k.analytic = delta_intra[tmp.idx2,i])    
  }
  
} else {
  
  sig_fm[] <- TRUE
  
  # evaluate at desired strikes 
  p_intra <- delta_intra <- gamma_intra <- matrix(NA,sum(internal1D),length(t_grid))
  for(j in 1:length(t_grid)){
    tmp  <- c.spline_predict(k_grid[internal1D], arbFreeIV$K[j,sig_fm[j,]], arbFreeIV$p[j,sig_fm[j,]], 
                             arbFreeIV$gamma[j,sig_fm[j,]])
    p_intra[,j] <- tmp$g
    delta_intra[,j] <- tmp$delta
    gamma_intra[,j] <- tmp$gamma
  }
  
  # Benaim EXTRAPOLATION:  
  # numerical derivatives 
  p_extra.numeric <- t(benaimEtrapol(t(p_intra),k_grid,internal1D,S0,r,q,t_grid))
  # analytical derivatives 
  p_extra <-  t(benaimEtrapol(t(p_intra),k_grid,internal1D,S0,r,q,t_grid, 
                              c_kk.analytic = gamma_intra, c_k.analytic = delta_intra))
  
}

# PLOT prices and partial derivatives 
tmp <- c(min(k_grid),2500)
matplot(p_extra, x=k_grid, type='l', main='extrapolated prices',xlim=tmp)
matplot(diffMat(K,1)$mat%*%t(p_mkt),x=K,type='l',main='Initial P_k')
matplot(diffMat(k_grid,1)$mat%*%p_extra.numeric,x=k_grid,type='l', 
        main='P_k (numerical)',xlim=tmp)
matplot(diffMat(k_grid,1)$mat%*%p_extra,x=k_grid,type='l',
        main='P_k (analytical)',xlim=tmp)
matplot(diffMat(K,2)$mat%*%t(p_mkt),x=K,type='l',main='Initial P_kk')
matplot(diffMat(k_grid,2)$mat%*%p_extra.numeric,x=k_grid,type='l', 
        main='P_kk (numerical)',xlim=tmp)
matplot(diffMat(k_grid,2)$mat%*%p_extra,x=k_grid,type='l', 
        main='P_kk (analytical)',xlim=tmp)

# compute IV 
iv_extra <- matrix(NA,nk,nt)
for(i in 1:nk){
  for(j in 1:nt){
    iv_extra[i,j] <- BScall_iv(p_extra[i,j],S0,k_grid[i],r[j],q[j],t_grid[j])
  }
}


# -- 1.2. Simple TPS ---------------------------------------------------------------------

IVbig <- matrix(NA,length(k_grid),length(t_grid))

# total variance intrapolation  
xy <- expand.grid(mat, K)
tpsFit <- Tps(xy, as.vector(mat*(iv^2)), 
              lambda=1e-10, ### <<<<---- set lambda 
              m=3)          ### <<<<---- set degree 
new.xy <- expand.grid(t_grid, k_grid[internal1D])
l <- lapply(1:length(t_grid), FUN=function(j){cbind(t_grid[j],k_grid[internal1D])})
new.xy <- do.call(rbind, l)
IVbig[internal2D] <- sqrt(predict(tpsFit,x=new.xy)/new.xy[,1])

# smooth, gradual flattening of the volatility curve
IVbig[,] <- t(flateningEtrapol(x=t(IVbig[internal1D,]), grid=log(k_grid), internal=internal1D, 
                               method='exponential', 
                               lambda=2)) ### <<<<---- set lambda 

# compute call prices 
p.directExtrap <- t(BScall_price2D(S0, k_grid, r, q, t(IVbig), t_grid)) 

# PLOT prices and partial derivatives 

tmp <- c(min(k_grid),2500)
matplot(p.directExtrap,x=k_grid,type='l',xlim=tmp)
matplot(diffMat(k_grid,1)$mat%*%p.directExtrap,x=k_grid,type='l',xlim=tmp)
matplot(diffMat(k_grid,2)$mat%*%p.directExtrap,x=k_grid,type='l',xlim=tmp)

# 2. Local Volatility ------------------------------------------------------------------------

# which LV surface? 
# p_extra <- p.directExtrap

######### LV boudreis 
max.lv <- 1.5
min.lv <- 0.04
nonSig.lv <- mean(data[[t]]$iv[,6])

# -- 2.1. ANDERSEN 97 ------------------------------------------------------------------------

abr <- AndersenBrothertonR(p_extra, s_grid=k_grid, t_grid, s.idx, r, q, theta=0.5, 
                           constrain=NULL, sig.idx=sig, 
                           imposeBoundryValues='grid_boundry',  
                           spotRates=T, ns.iv=nonSig.lv, max.iv=max.lv)

lv_andersen <- sqrt(abr$a*2)

open3d()  
visualizeSurface(log(k_grid), c(t_grid), lv_andersen )
open3d()
visualizeSurface(log(k_grid), c(0,t_grid), abr$A_ini )

# -- 2.2. DUPIRE 97 (price formula) ---------------------------------------------------------  

rf <- spotTOforward(r,t_grid,'contTOcont')$f
qf <- spotTOforward(q,t_grid,'contTOcont')$f

c_k <- diffMat(k_grid,1)$mat%*%p_extra
c_kk <- diffMat(k_grid,2)$mat%*%p_extra
c_t <- t(diffMat(t_grid,1)$mat%*%t(p_extra))
c_kk.ini <- (diffMat(K,2)$mat%*%t(p_mkt))

lv.dup <- t(dupireLV_p(t(p_extra),k_grid,rf,qf,t(c_t),t(c_k),t(c_kk)))

lv.dup[!sig] <- nonSig.lv
lv.dup[lv.dup>max.lv] <- max.lv
lv.dup[lv.dup<min.lv | is.na(lv.dup)] <- min.lv

# -- 2.3. DUPIRE 97 (IV formula) ------------------------------------------------------------

sig_k <- diffMat(k_grid,1)$mat%*%IVbig
sig_kk <- diffMat(k_grid,2)$mat%*%IVbig
sig_t <- t(diffMat(t_grid,1)$mat%*%t(IVbig))

lv.dup2 <- t(dupireLV_iv(t(IVbig),k_grid,rf,qf,t(sig_t),t(sig_k),t(sig_kk),t_grid,S0))

lv.dup2[!sig] <- nonSig.lv
lv.dup2[lv.dup2>max.lv] <- max.lv
lv.dup2[lv.dup2<min.lv | is.na(lv.dup2)] <- min.lv

# PLOTS:
# compare Dupire LVS compute with the two formulas 
open3d()  
layout3d(matrix(1:2,1,2) , sharedMouse = TRUE)
visualizeSurface(log(k_grid), t_grid, (lv.dup2),title='Dupire (IV formula)')
next3d()  
visualizeSurface(log(k_grid), t_grid, (lv.dup),title='Dupire (price formula)')

# Plot Dupire vs Andersen 
open3d()  
layout3d(matrix(1:2,1,2) , sharedMouse = TRUE)
visualizeSurface(log(k_grid), t_grid, (lv.dup2),title='Dupire')
next3d()  
visualizeSurface(log(k_grid), t_grid, (lv_andersen),title = 'Andersen')


# 3. PRICING ----------------------------------------------------------------------------

andersen.fdm.p <- rep(NA,length(K))
for(i in 1:length(K)){
  andersen.fdm.p[i] <- fdmABR(K[i], s_grid=k_grid, s.idx=s.idx, t_grid, abc=abr, r, q, 
                              theta = 0.5, n.imp = 0, spotRates = TRUE)$p
}
dupire.fdm.p <- rep(NA,length(K))
for(i in 1:length(K)){
  dupire.fdm.p[i] <- fdmBack(K[i], k_grid, s.idx, t_grid, lv.dup, r, q, theta = 0.5, 
                             n.imp = 1, ab.adjust = TRUE, spotRates = TRUE)$p
}
dupire2.fdm.p <- rep(NA,length(K))
for(i in 1:length(K)){
  dupire2.fdm.p[i] <- fdmBack(K[i], k_grid, s.idx, t_grid, lv.dup2, r, q, theta = 0.5, 
                             n.imp = 1, ab.adjust = TRUE, spotRates = TRUE)$p
}

cbind(
  K = K/S0,
  mkt_p = p_mkt[mat.idx,],
  # Andersen 
  And_p = andersen.fdm.p,
  diff = round((andersen.fdm.p - p_mkt[mat.idx,]),2),
  relDiff = round((andersen.fdm.p/p_mkt[mat.idx,]-1)*100,2),
  # Dupire Price Formula 
  dup.p_p = dupire.fdm.p,
  diff = round((dupire.fdm.p - p_mkt[mat.idx,]),2),
  relDiff = round((dupire.fdm.p/p_mkt[mat.idx,]-1)*100,2),
  # Dupire IV Formula
  Dup.iv_p = dupire2.fdm.p,
  diff = round((dupire2.fdm.p - p_mkt[mat.idx,]),2),
  relDiff = round((dupire2.fdm.p/p_mkt[mat.idx,]-1)*100,2)
)
print('RMSE:')
cbind(
  sqrt(mean((p_mkt[mat.idx,] - andersen.fdm.p)^2)),
  sqrt(mean((p_mkt[mat.idx,] - dupire.fdm.p)^2)),
  sqrt(mean((p_mkt[mat.idx,] - dupire2.fdm.p)^2))  
)

