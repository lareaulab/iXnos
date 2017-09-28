# leave-one-out analysis of codon position contributions
# supp figure B
# data: Green yeast

args <- commandArgs(trailingOnly = TRUE)
mses_fname = args[1]
out_fname = args[2]

dt = read.delim(mses_fname, header = T, skip = 1)

leaveout.mean = apply(dt[1:13,2:11], 1, mean)
full.mean = apply(dt[14,2:11], 1, mean)

leaveout.se = apply(dt[1:13,2:11], 1, sd) / sqrt(10)
full.se = apply(dt[14,2:11], 1, sd) / sqrt(10)

diff = leaveout.mean - full.mean
diff.se = sqrt( leaveout.se^2 + full.se^2 )

label1 = c(-7,NA,-5,NA,-3,NA,"P",NA,1,NA,3,NA,5)
label2 = c(NA,-6,NA,-4,NA,"E",NA,"A",NA,2,NA,4,NA)

ymin = min( diff - 2 * diff.se )
ymax = max( diff + 2 * diff.se )

pdf( out_fname, width=2, height=1.67, pointsize = 7, useDingbats = F, bg = "white" )
#cairo_pdf( out_fname, width=2, height=1.67, pointsize = 7)
par( mex = 0.65 )
par( mar = c(6,5.5,5,3) )
par( oma = c(0,1.5,1,0) )
par( lwd = 0.75 )
centers = barplot( diff,
                   ylim = c( ymin, max(ymax, round(ymax,1)) ),
                   col = "mediumpurple3",
                   border = "mediumpurple3",
                   space = 0.3,
                   axes = F,
                   axisnames = F,
                   ylab = "MSE increase",
                   xlab="codon position"
)
axis( 2, at = seq( round( ymin, 1 ), round( ymax, 1 ), by = 0.1 ), lwd = 0.75 )
axis( 1, at = centers, padj = -1,
      labels =label1,
      tick=F, cex.axis = 0.7)
axis( 1, at = centers, padj = -1,
      labels = label2,
      tick=F, cex.axis = 0.7)
segments( centers, diff - diff.se, 
          centers, diff + diff.se )
arrows( centers, diff - diff.se, 
        centers, diff + diff.se,
        angle = 90, code = 3, 
        length = 0.02 )
mtext( "B", font = 2, line = -3, side = 3, outer = T, adj = 0 )
dev.off()

