# binned error of true and predicted scaled counts
# figure 1E
# data: Weinberg yeast

args <- commandArgs(trailingOnly = TRUE)
y_te_fname = args[1] 
y_te_hat_fname = args[2]
out_fname = args[3]

y_te = read.delim(y_te_fname, header=F)
y_te_hat = read.delim(y_te_hat_fname, header=F)

pts = length(y_te$V1)
bins = 700
binsize = round(pts/bins)

true_count_order = order(y_te$V1)

y_te_ordered = y_te$V1[true_count_order]
y_te_hat_ordered = y_te_hat$V1[true_count_order]

sq_errs = (y_te_ordered - y_te_hat_ordered)^2

avg_errs = sapply( seq(1, pts, by = binsize), function(x){sum(sq_errs[x:(x+binsize)])/binsize })

# scaled count value at the boundary of each bin
bin_cutoffs = y_te_ordered[ seq(1, pts, by = binsize) ]

cairo_pdf( out_fname, width=2, height=1.67, pointsize = 7 )
par( mex = 0.65 )
par( mar = c(6,4.5,5,3) )
par( oma = c(0,1.5,1,0))
par( lwd = 0.75 )

plot(1:length(avg_errs), log(avg_errs)/log(10), 
     pch = 20, cex = 0.5, 
     bty = "n", 
     axes = F,
     ylim = c(-1,2),
     xlab = "bin number",
     ylab = "log MSE",
     col="cornflowerblue")

axis(1, at = seq(0, length(bin_cutoffs), by = 100),
     labels = seq(0, length(bin_cutoffs), by = 100),
     lwd = 0.75)
axis(2, at = seq(-1,2,by=1), lwd = 0.75)
axis(3, at = seq(0, length(bin_cutoffs), by = 100), 
     labels = round(bin_cutoffs[seq(1, length(bin_cutoffs), by = 100)], 1),
     lwd = 0, lwd.ticks = 0.75, cex.axis = 0.7, padj = 0.5
)
mtext( "scaled count", side = 3, line = 2, cex = 0.7)
mtext( "F", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 

dev.off()