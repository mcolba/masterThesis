
###  Fengler_example.R  ### 

rm(list = ls())

source("./source/BSformulas.R")
source("./source/plotter.R")
source("./source/LV/pre_smoother.R")
source("./source/LV/Fengler09.R")
source("./source/LV/cubicSpline.R")
source("./source/LV/Dupire94.R")


# load dataset  
load("./data/data.RData")


# EXAMPLE 1.1 ---------------------------------------------------------------------------
# Pre-Smoothing: old vs new grid 

#  set date 
t=1
print(data[[t]]$date)
# initial expiration and forward-moneyness coordinates 
tmp <- expand.grid(data[[t]]$forward, data[[t]]$strikes)
moneyness <- tmp[,2]/tmp[,1]
initialXY <- cbind.data.frame(expiry = expand.grid(data[[t]]$expiries, data[[t]]$strikes)[,1], 
                              moneyness)

# set forward-moneyness grid 
nk <- 10 
k_grid <- seq(0.75, 1.2, length=nk)
add_t <- 14/12
nt <- length(c(data[[t]]$expiries, add_t))


# pre-smoothe the ivs 
arguments <- c(data[[t]], list(k_grid = k_grid,add_t = add_t))
preSmooth <- do.call(TPSpreSmoother, arguments)

# extract the new (regular) grid in foward moneynes 
newGrid <- preSmooth$xy
# grups for ggplot 
newGrid$g1 <- rep(1:nt, nk) 
newGrid$g2 <- sort(rep(1:nk, nt))

# 2D plot -------------------------------------------------------------------------------

library(ggplot2)
ggplot(newGrid, aes(x=expiry, y=k)) +  
  geom_point(alpha = 0.5) + 
  geom_line(aes(group = g1), linetype="dotted") +
  geom_line(aes(group = g2), linetype="dotted") +
  geom_point(data = initialXY, aes(x=expiry, y=moneyness), col = "red", alpha = 0.5)

# PLOT presmoothed IVS (forward moneyness) ----------------------------------------------
xyz <- cbind(initialXY, iv = as.vector(data[[t]]$iv))
open3d()
points3d(xyz, size = 7) # old knots 
visualizeSurface(preSmooth$t_grid , preSmooth$k_grid, preSmooth$iv,
                 zlab='IV',xlab='time',ylab='K/S') 

# PLOT presmoothed price surface (forward moneyness) ------------------------------------
initial_p <- as.vector(do.call(BScall_price2D, data[[t]]))
xyz <- cbind(initialXY, price = initial_p)
open3d()
points3d(xyz, size = 7) # new knots 
visualizeSurface(preSmooth$t_grid , preSmooth$k_grid, preSmooth$p) # original surface 


# EXAMPLE 1.2 ---------------------------------------------------------------------------
# quadratic program 

lambda = 1e-2
# solve quadratic programm 
arguments2 <- c(preSmooth, list(smooth = lambda, spot = data[[t]]$spot))
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
points(preSmooth$p[exp,], x = cs$K[exp,], type = "p") # from smoothed greed 
initial_p <- do.call(BScall_price2D, data[[t]])
points(initial_p[exp,], x = data[[t]]$strikes, type = "p", pch = 8) # from smoothed greed 

# PLOT Arbitrage-Free price surface -----------------------------------------------------
open3d()
visualizeSurface(x = preSmooth$t_grid, y = points, z = val, alpha=0.3) # interpolated points 
for (i in 1:nrow(cs$p)){
  lines3d(x = preSmooth$t_grid[i], y = points, z = val[i,], lwd = 3, col = "green")
}
xy <- expand.grid(data[[t]]$expiries, data[[t]]$strikes)
xyz <- cbind(xy, as.vector(initial_p))
points3d(xyz, size = 6) # original prices 




