
###  examples.R  ### 

rm(list = ls())

source("./source/BSformulas.R")
source("./source/plotter.R")
source("./source/LV/pre_smoother.R")
source("./source/LV/Fengler09.R")
source("./source/LV/cubicSpline.R")
source("./source/LV/Dupire.R")


# load dataset  
load("./Rdata/data.RData")

#
# EXAMPLE 1
# Pre-Smoothing: old vs new grid 
#

#  set date 
t=1
# initial expiration and forward-moneyness coordinates 
tmp <- expand.grid(data[[t]]$forward, data[[t]]$strikes)
moneyness <- tmp[,2]/tmp[,1]
initialXY <- cbind.data.frame(expiry = expand.grid(data[[t]]$expiries, data[[t]]$strikes)[,1], 
                              moneyness)
# set forward-moneyness grid 
nk <- 10 
k_grid <- seq(0.75, 1.2, length=nk)
add_t <- c(0.8, 0.9)
nt <- length(c(data[[t]]$expiries, add_t))

# pre-smoothe the ivs 
arguments <- c(data[[t]], list(k_grid = k_grid,
                               add_t = add_t))
preSmooth <- do.call(TPSpreSmoother, arguments)

# extract the new (regular) grid in foward moneynes 
newGrid <- preSmooth$xy
# grups for ggplot 
newGrid$g1 <- rep(1:nt, nk) 
newGrid$g2 <- sort(rep(1:nk, nt))

# 2D plot 
library(ggplot2)
ggplot(newGrid, aes(x=expiry, y=k)) +  
  geom_point(alpha = 0.5) + 
  geom_line(aes(group = g1), linetype="dotted") +
  geom_line(aes(group = g2), linetype="dotted") +
  geom_point(data = initialXY, aes(x=expiry, y=moneyness), col = "red", alpha = 0.5)

# 3D IVS plot (forward moneyness)
xyz <- cbind(initialXY, iv = as.vector(data[[t]]$iv))
open3d()
points3d(xyz, size = 7) # old knots 
visualizeSurface(preSmooth$t_grid , preSmooth$k_grid, preSmooth$iv) # new surface 

# 3D prices plot (forward moneyness)
initial_p <- as.vector(do.call(BScall_price2D, data[[t]]))
xyz <- cbind(initialXY, price = initial_p)
open3d()
points3d(xyz, size = 7) # new knots 
visualizeSurface(preSmooth$t_grid , preSmooth$k_grid, preSmooth$p) # original surface 

#
# EXAMPLE 2 
# quadratic program
#

# rm(list = ls())

# set parameters 
t=1
k_grid_length = 50
lambda = 1e-2
add_t = c(0.79, 0.83, 0.87, 0.90, 0.93, 0.96, 0.99, 1)

# presmooth the ivs 
k_grid <- seq(0.70, 1.25, length = k_grid_length)
arguments <- c(data[[t]], list(k_grid = k_grid,
                               add_t = add_t))
preSmooth <- do.call(TPSpreSmoother, arguments)

# solve quadratic programm 
arguments2 <- c(preSmooth, list(smooth = lambda,
                                spot = data[[t]]$spot))
cs <-do.call(solveQuadprog, arguments2) 

# compute intermediate strikes 
points <- seq(min(data[[t]]$strikes)-0, max(data[[t]]$strikes)+0, length.out = 100)
val <- matrix(NA, length(preSmooth$t_grid), 100)
for (i in 1:100){
  for (j in 1:length(preSmooth$t_grid)){
    val[j,i] <- evalSplin(points[i], cs$K[j,], cs$p[j,], cs$gamma[j,])
  }
}

# 2D 
exp = 1
plot(val[exp,], x = points, type = "l") # interpolated prices 
points(preSmooth$p[exp,(ncol(preSmooth$p):1)], x = cs$K[exp,], type = "p") # from smoothed greed 
initial_p <- do.call(BScall_price2D, data[[t]])
points(initial_p[exp,], x = data[[t]]$strikes, type = "p", pch = 8) # from smoothed greed 

# 3D 
open3d()
visualizeSurface(x = preSmooth$t_grid, y = points, z = val, alpha=0.3) # interpolated points 
for (i in 1:nrow(cs$p)){
  lines3d(x = preSmooth$t_grid[i], y = points, z = val[i,], lwd = 3, col = "green")
}
xy <- expand.grid(data[[t]]$expiries, data[[t]]$strikes)
xyz <- cbind(xy, as.vector(initial_p))
points3d(xyz, size = 6) # original prices 


## NEW GRID ##
K2 <- seq(data[[t]]$spot*0.8, data[[t]]$spot*1.2, length.out = k_grid_length)
cs2 <- c.spline_predict(K2 ,cs$K ,cs$p, cs$gamma)

# LOCAL VOL SURFACE 
# to improve ?
nt <- length(preSmooth$t_grid)
c_t <- (cs2$g[2:nt,] - cs2$g[1:(nt-1),])/
  (preSmooth$t_grid[2:nt]-preSmooth$t_grid[1:(nt-1)])
# Reduce grid
gamma <- cs2$gamma[1:(nt-1),]
p <- cs2$g[1:(nt-1),]
delta <- cs2$delta[1:(nt-1),]
t_grid <- preSmooth$t_grid[1:(nt-1)]
rates <- preSmooth$rates[1:(nt-1)]
dividends <- preSmooth$dividends[1:(nt-1)]

LV <- dupireLV_p(p, K2, c_t, delta, gamma, rates, dividends)
LV[LV>1]<- 1

idx_na <- which(is.na(LV), arr.ind =  T)
# LV[, c(1,ncol(LV))] <- LV[, c(2, ncol(LV)-1)]

idx <- which(is.na(LV), arr.ind =  T)
# idx <- matrix(idx[idx[,2]!=1 & idx[,2]!=ncol(LV),], ncol=2)
if (ncol(idx)>1){
  u <- sort(unique(idx[,1]), decreasing = T)
} else u <- idx[,1]
for (i in u){
  idxk <- idx[idx[,1]==i, 2]
  LV[i, idxk] <- LV[(i+1), idxk]
}
xyz <- cbind(t_grid[idx[,1]], K2[idx[,2]], as.vector(LV[idx])) 

open3d()
visualizeSurface(x = t_grid, y = K2, z = LV)
points3d(xyz, size = 10, col = "red")


# Goodness of fit -> OK 
grid_p <- c.spline_predict(data[[t]]$strikes, cs$K, cs$p, cs$gamma)$g # [-c(6,7),]
(initial_p - grid_p)*100 # in cents 
(initial_p - grid_p)/grid_p * 100 # percent 


