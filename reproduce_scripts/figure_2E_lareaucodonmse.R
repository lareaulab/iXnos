# leave-one-out analysis of codon position contributions
# figure 2E
# data: Lareau / Graham yeast

args <- commandArgs(trailingOnly = TRUE)
mses_fname = args[1]
out_fname = args[2]

dt = read.delim(mses_fname, header = T, comment = "#")

leaveout.mean = apply(dt[1:13,2:11], 1, mean)
full.mean = apply(dt[14,2:11], 1, mean)

data = data.frame( cbind( t(dt[14,2:11]), NA, t(dt[1:13,2:11]) ))

label1 = c("full", NA, -7, NA, -5, NA, -3, NA,  "P", NA,  1,  NA, 3,  NA, 5)
label2 = c(NA,     NA, NA, -6, NA, -4, NA, "E", NA,  "A", NA, 2,  NA, 4,  NA)

ymin = floor( min(data, na.rm=T) * 20 ) / 20
ymax = ceiling( max(data, na.rm=T) * 20 ) / 20

pdf( out_fname, width=2, height=1.67, pointsize = 7, useDingbats = F, bg = "white" )
#cairo_pdf( out_fname, width=2, height=1.67, pointsize = 7 )
par( mex = 0.65 )
par( mar = c(6,5.5,4,3) )
par( oma = c(0,1.5,1,0) )
par( lwd = 0.75 )

plot( NA, 
      xlim = c( 1, 15 ),
      ylim = c( ymin, ymax ),
      axes = F,
      xlab = "codon position",
      ylab = "MSE"
)

rect( c(2.6:14.6), rep(full.mean,13), c(3.4:15.4), leaveout.mean, 
      col = "mediumpurple3", border = NA )

stripchart( data, 
            vertical = T, 
            pch = 16,
            col = rgb(0.2,0.2,0.2,.3),
            cex = .4,
            add = T
            )

axis( 2, lwd = 0.75 )
axis( 1, at = 1:15, padj = -1, labels = label1, tick = F, cex.axis = 0.7)
axis( 1, at = 1:15, padj = -1, labels = label2, tick = F, cex.axis = 0.7)

abline( h = full.mean, lty = 3, lwd = 0.5 )

mtext( "E", font = 2, line = -3, side = 3, outer = T, adj = 0 )
dev.off()

