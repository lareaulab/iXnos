# zoomed in scatterplot of true and predicted scaled counts
# figure 1F
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

pdf( out_fname, width=2, height=1.67, pointsize = 7, useDingbats = F, bg = "white" )
par( mex = 0.65 )
par( mar =c(5, 5.5, 2, 3))
par( oma = c(0,1.5,1,0) )
par( lwd = 0.75 )
plot( NA, 
      xlim = c( 0, 2 ),
      ylim = c( 0, 2 ),
      axes = F,
      xaxs = "i",
      bty = "n",
      asp = 1,
      xlab = "true scaled count",
      ylab = "predicted scaled count")
smoothScatter( y_te$V1, 
               y_te_hat$V1, 
               colramp = jet.colors.alpha, 
               nbin=c(1024, round(1024*aspratio)), 
               bandwidth = 0.1, 
               add = T,
               pch=18, cex=0.4, col="darkgray",
               nrpoints = 100,
               postPlotHook = NULL)
axis( 1, at = 0:2, lwd = 0.75 )
axis( 2, at = 0:2, lwd = 0.75 )
mtext( "F", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 
dev.off()