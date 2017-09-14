# scatterplot of true and predicted scaled counts
# figure 1D
# data: Weinberg yeast

library("fields")

options("preferRaster" = T) # the scale bar and plot will both be raster images 

args <- commandArgs(trailingOnly = TRUE)
y_te_fname = args[1] 
y_te_hat_fname = args[2]
out_fname = args[3]

y_te = read.delim(y_te_fname, header=F)
y_te_hat = read.delim(y_te_hat_fname, header=F)

jet.colors.alpha <- colorRampPalette(c("#00007F00",
                                        "#0000FF00", 
                                        "#007FFF30", 
                                        "#00FFFF50",
                                        "#7FFF7F80", 
                                        "#FFFF00a0", 
                                        "#FF7F00c0", 
                                        "#FF0000e0", 
                                        "#AA0000ff"
                                        ), 
                                     alpha=T,
                                     interpolate = "spline"
                                     )

aspratio = max(y_te_hat) / max(y_te)

scalebar <- function(){
  zvals <- get('dens', envir = parent.frame(1))
  zlim <- c( min(zvals, na.rm=T), max(zvals, na.rm=T) )
  colramp <- get('colramp', parent.frame(1))
  fields::image.plot( zlim = zlim, col = colramp(256), 
                      legend.only = T, add = F, 
                      axis.args = list( lwd = 0.75))
}


cairo_pdf( out_fname, width=6, height=1 + 5*aspratio, pointsize=7 )
par( mex = 0.65 )
par( mar =c(6,5.5,5,7) )
par( oma = c(0,1.5,1,0 ))
par( xpd = T )
par( lwd = 0.75 )
# a hack to add the axes I want because of a bug in image.plot that defines a
# truncated plot window
plot( NA, xlim = c( 0, max(y_te$V1) ), ylim = c( 0, max(y_te_hat$V1) ),
      axes = F,
      bty = "n",
      asp = 1,
      xlab = "true scaled count",
      ylab = "predicted scaled count"
      )
axis( 1, pos = -1, at = c(0:5,10,20,30), lwd = 0.75 )
axis( 2, pos = -1, at = 0:5, lwd = 0.75 )

mtext( "D", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 

smoothScatter( y_te$V1, y_te_hat$V1, 
               colramp=jet.colors.alpha, 
               nbin=c(512, round(512*aspratio)), 
               bandwidth = 0.1, 
               add = T,
               pch=18, cex=0.4, col="darkgray",
               nrpoints=100, 
               postPlotHook=scalebar )


dev.off()
