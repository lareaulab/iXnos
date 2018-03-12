# downsampling to different coverage
# supp figure, panel A
# full model on all genes

args <- commandArgs(trailingOnly = TRUE)
corrs_fname = args[1] # corrs_by_gene_density.txt
out_fname = args[2]

corrs = read.delim( corrs_fname, header=T, comment.char = "#", row.names=1 )

num = length( row.names(corrs) )
training = which( as.logical(corrs$tr_gene) )
nottraining = which( !as.logical(corrs$tr_gene) )

labels = c("all data",
           "subsampled to match\n1000th gene",
           "subsampled to match\n2000th gene",
           "subsampled to match\n3000th gene",
           "subsampled to match\n4000th gene"
           )

plot_corrs = function( x, y, ylim1, ylim2 ) { 
  plot( 1:num,
        x,
        axes = F,
        pch = NA, 
        bty = "n",
        ylim = c(ylim1, ylim2),
        xlab = NA,
        ylab = NA
  )
  points( training, x[training], bg = "#0000ff25", col = NA, pch = 21, cex = 0.5 )
  points( nottraining, x[nottraining], bg = "#00000040", col = NA, pch = 21, cex = 0.5 )
  axis(2, at = c(0,0.5,1), labels = c(0, NA, 1), las = 1)
  loess.line = predict( loess( x[nottraining] ~ nottraining ),
                        1:num )
  lines( 1:num, loess.line, col="red", lwd = 1.5 )
  mtext( y, side = 3, line = -2.5, cex = 0.75, adj = 1 )
}


pdf( out_fname, width=2, height=3.5, pointsize=7, useDingbats=F)
par( mfrow = c(5,1) )
par( mex = 0.65 ) # sets margin stuff
par( lwd = 0.75 )
par( xpd = NA )
par(oma = c(5, 1.5, 6, 0),
    mar = c(1, 5.5, 0, 3) )
par( cex = 1 )

Map( plot_corrs, corrs[,4:8], labels, 0, 1 )
axis( 1, at = seq(0, 4000, by=1000), labels = c("0",NA,"2000",NA,"4000") )
#axis( 2, at = c(0,0.5,1) )
title( xlab = "genes sorted by footprint density" )
mtext( "Pearson correlation per gene", side = 2, outer = T, line = -2.5)
mtext( "A", font = 2, line = 2, side = 3, outer = T, adj = 0 ) 

dev.off()
