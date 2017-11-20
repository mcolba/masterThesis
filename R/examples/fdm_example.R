
### example_fdm.R ### 

rm(list = ls())

source("./source/numerical/fdm.R")
source("./source/LV/mesh.R")
source("./source/BSformulas.R")


#### EXAMPLE 1: strike independent IVS - pricing calls using Backward-FDM ####

# FINDINGS: 
# - the ab-adjustment improve accuracy in the Crank-Nicolson scheme only
# - non-equispaced mesh improve accuracy, use: alpha = 0.1 or lower
# - small x-steps near ATM produce disturbance in delta and prices -> 
#   use a mixed scheme or a bigger alpha    

# Call-Option parameters 
mat <- 2; S0 <- 590
K <- S0*seq(0.75,1.25,by=0.05) 

# FDM parameters
nk <- 150
nt <- 75
t_grid = seq(0.15,mat,length.out=nt)
dt <- c(t_grid[1], t_grid[2:nt]-t_grid[1:nt-1])

# mkt parameters  
r <- seq(0.020, 0.100, length.out=nt)
q <- seq(0.046, 0.016, length.out=nt)

# Strike independent implied volatility 
LV <- matrix(rep(seq(0.2, 0.1, length.out=nt), nk), nk, nt, byrow=T)
IV <- sqrt(cumsum(LV[1,]^2*dt)/t_grid)

# k-grid
min.s <- S0*exp((r[nt]-q[nt]-0.5*0.2^2)*mat - 5*0.2*sqrt(mat)) 
max.s <- S0*exp((r[nt]-q[nt]-0.5*0.2^2)*mat + 5*0.2*sqrt(mat)) 

# construct space-grid: 
# 1. Equally spaced #
k_grid <- seq(min.s, max.s, length.out=nk-1)
# nclude spot in the mesh 
k_grid <- c(k_grid, S0)
s_idx <- which(order(k_grid)==length(k_grid))
k_grid <- sort(k_grid)
# 2. hyperbolic # 
hyp <- hyperbolicMesh(min.s, max.s, c=S0, n=nk, alpha=0.001, x_star=S0)
k_grid <- hyp$x
s_idx <- hyp$idx

# analytical prices 
IV_2d <- matrix(rep(IV,length(K)),length(K),nt, byrow=T)
p_BS <- BScall_price2D(spot=S0, strikes=K, rates=r, dividends=q, iv=t(IV_2d), expiries=t_grid)

# declare 
back_cn_dicrete <- back_mix_dicrete <- back_imp_dicrete <- back_cn_cont <- 
  back_mix_cont <- back_imp_cont <- matrix(NA,length(K))

# Backward FDM 
for(i in 1:length(K)){
  
  back_cn_dicrete[i] <- fdmBack(K=K[i], s_grid=k_grid, s.idx=s_idx, t_grid=t_grid, lv=LV, r=r, 
                                q=q, theta = 0.5, n.imp = 0, ab.adjust = TRUE, spotRates = TRUE)$p 
  back_mix_dicrete[i] <- fdmBack(K=K[i], s_grid=k_grid, s.idx=s_idx, t_grid=t_grid, lv=LV, r=r, 
                                 q=q, theta = 0.5, n.imp = 2, ab.adjust = TRUE, spotRates = TRUE)$p 
  back_imp_dicrete[i] <- fdmBack(K=K[i], s_grid=k_grid, s.idx=s_idx, t_grid=t_grid, lv=LV, r=r, 
                                 q=q, theta = 0, n.imp = 0, ab.adjust = TRUE, spotRates = TRUE)$p 
  back_cn_cont[i] <- fdmBack(K=K[i], s_grid=k_grid, s.idx=s_idx, t_grid=t_grid, lv=LV, r=r, 
                             q=q, theta = 0.5, n.imp = 0, ab.adjust = FALSE, spotRates = TRUE)$p 
  back_mix_cont[i] <- fdmBack(K=K[i], s_grid=k_grid, s.idx=s_idx, t_grid=t_grid, lv=LV, r=r, 
                              q=q, theta = 0.5, n.imp = 2, ab.adjust = FALSE, spotRates = TRUE)$p 
  back_imp_cont[i] <- fdmBack(K=K[i], s_grid=k_grid, s.idx=s_idx, t_grid=t_grid, lv=LV, r=r, 
                              q=q, theta = 0, n.imp = 0, ab.adjust = FALSE, spotRates = TRUE)$p 
}

table <- cbind(
  c( p_BS[nt,], 0 ),
  c( back_cn_dicrete-p_BS[nt,], mean((back_cn_dicrete-p_BS[nt,])^2) ),
  c( back_mix_dicrete-p_BS[nt,], mean((back_mix_dicrete-p_BS[nt,])^2) ),
  c( back_imp_dicrete-p_BS[nt,], mean((back_imp_dicrete-p_BS[nt,])^2) ),
  c( back_cn_cont-p_BS[nt,], mean((back_cn_cont-p_BS[nt,])^2) ),
  c( back_mix_cont-p_BS[nt,], mean((back_mix_cont-p_BS[nt,])^2) ),
  c( back_imp_cont-p_BS[nt,], mean((back_imp_cont-p_BS[nt,])^2) )
)
colnames(table) <- c('BS','bk.CN.disc','bk.MIX.disc','bk.IMP.disc','bk.CN.cont','bk.MIX.cont','bk.IMP.cont')
rownames(table) <- c(round(K/S0,2), 'MSE') 
print(table)

#### EXAMPLE 2: strike independent IVS, ATM delta with Backward-FDM ####

cn <- fdmBack(K=K[6], s_grid=k_grid, s.idx=s_idx, t_grid=t_grid, lv=LV, r=r, 
                              q=q, theta = 0.5, n.imp = 0, ab.adjust = TRUE, spotRates = TRUE) 
mix <- fdmBack(K=K[6], s_grid=k_grid, s.idx=s_idx, t_grid=t_grid, lv=LV, r=r, 
                              q=q, theta = 0.5, n.imp = 1, ab.adjust = TRUE, spotRates = TRUE) 
imp <- fdmBack(K=K[6], s_grid=k_grid, s.idx=s_idx, t_grid=t_grid, lv=LV, r=r, 
                              q=q, theta = 0, n.imp = 0, ab.adjust = TRUE, spotRates = TRUE)
# plot 
idx <- cn$s_grid>300 & cn$s_grid<800
plot((diffMat(k_grid)$mat%*%cn$v)[idx], x=cn$s_grid[idx], type="l", lty=6, ylab = "Delta", xlab = "Spot")
lines((diffMat(k_grid)$mat%*%imp$v)[idx], x=imp$s_grid[idx], col="green", lty = 5)
lines((diffMat(k_grid)$mat%*%mix$v)[idx], x=mix$s_grid[idx], lty=2, col = "red")
legend("bottomright",legend=c("Crank-Nicolson","Implicit","Mixed"), 
       lty=c(6,5,2), col = c("black","green","red"))

