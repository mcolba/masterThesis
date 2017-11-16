
###  plotter.R  ###

# load/install packages 
if (!require("rgl")) install.packages("rgl") # for 3D plots

# matlab-like color scale  
matColors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", 
                                    "yellow", "#FF7F00", "red", "#7F0000"))
g_colorTable <- matColors(300)  # set number of different colours
rm(matColors) 

# My default view point
g_myView <- matrix( 
  c(0.7448956, 0.6671253, -0.008632005, 0,
    -0.1403858, 0.1693734, 0.975502133, 0,
    0.6522444, -0.7254353, 0.219820470, 0,
    0, 0, 0, 1), byrow = TRUE, ncol = 4 )

# zoom<-par3d()$zoom
# userMatrix<-par3d()$userMatrix
# windowRect<-par3d()$windowRect

visualizeSurface_iv <- function(day){
  #
  # Visualize surface given date
  #
  
  i <- day 
  z <- data[[i]]$iv 
  x <- as.numeric(data[[i]]$expiries)
  y <- as.numeric(log(data[[i]]$spot/data[[i]]$strikes)) # ln(S/K)
  # call to the general visualizeSurface functon 
  return(visualizeSurface(x, y, z, title = as.character(data[[i]]$date), 
                          zlab = 'IV', xlab = 'Expiration', ylab = 'ln(S/K)',
                          c.min = 0, c.max = 0.9))
}

visualizeSurface <- function(x, y, z, title = '', zlab = '', xlab = '', ylab = '', zm=1.15, 
                             c.min = min(z, na.rm = T), c.max = max(z, na.rm = TRUE), ...){
  
  # check dimension 
  if(length(x)!=nrow(z)|length(y)!=ncol(z)){stop('wrong dimension')}
  
  arg <- list(...)
  
  # assign colours to different values of z 
  col <- g_colorTable[ findInterval(z, seq(c.min, c.max, length=300))]
  
  surface3d(x, y, z, color = col, smooth = FALSE, arg)
  aspect3d(1,1,1) # make it a square 
  axes3d(edges=c("x--", "y+", "z"), cex = 0.75) # axes position 
  title3d(main = title, cex = 1.1) # title 
  title3d(xlab = xlab, zlab = zlab, line = 2.3, cex = 0.9)  
  mtext3d(text = ylab, edge = "y+", line = 2.3, cex = 0.9) 
  grid3d(side = c("x","y+","z")) # grids position  

  # use my default definition, zoom, and view point
  par3d(zoom = zm, userMatrix = g_myView)
  par3d(windowRect=c(0,30,795,825))
  
}

legendBar <- function(min, max, nticks = 10, bottom = 1, left = 1,
                      top = bottom, right = left, title = '',
                      colorTable = g_colorTable, horizontal = TRUE){
  #
  #  plot a coloured bar (for the legend) 
  #  
  ticks <- round(seq(min, max, len=nticks),1) 
  scale <- (length(colorTable)-1)/(max-min)
  # dev.new(width=10, height=10)
  par(mar=c(bottom, left, top, right))
  if (horizontal==TRUE){
    plot(c(min,max), c(0, 10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
    for (i in 1:(length(colorTable)-1)) {
      y = (i-1)/scale + min
      rect(y,0,y+1/scale,10, col=colorTable[i], border=NA)
    }
  } else {
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
    for (i in 1:(length(colorTable)-1)) {
      y = (i-1)/scale + min
      rect(0,y,10,y+1/scale, col=colorTable[i], border=NA)
    }
  }
  axis(1, ticks, las=1)
  title(main = title, line = 1, font.main = 1, cex.main = 1.1)
}

saveCurretSurface <- function(fileName, square = FALSE){
  #
  # Save Surface in .png format 
  # 
  if (square == TRUE) par3d(windowRect=c(0,0,850,850))  # max definition on my screen 
  
  rgl.snapshot( paste("../../Editing/figures/", fileName, ".png", sep = ""), 
                fmt="png", top=TRUE)
}
