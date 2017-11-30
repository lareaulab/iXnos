# MSE by model type and sequence neighborhood 
# figure 1C
# data: Weinberg yeast

args <- commandArgs(trailingOnly = TRUE)
nn_mses_fname = args[1] #feat_neighborhood_mses.txt
lr_mses_fname = args[2] #linreg_mses.txt
str_mses_fname = args[3] #struc_mses.txt
out_fname = args[4]

full = read.delim(nn_mses_fname, header = T, comment = "#", row.names = 1)
lr = read.delim(lr_mses_fname, header = T, row.names = 1)
struc = read.delim(str_mses_fname, header = T, row.names = 1)

series_names = c("NN: codons", "NN: codons + nt", "LR: codons", "LR: codons + nt")
nn_cod = full$mean_MSE[c(1,3,5,7)]
nn_cod_nt = full$mean_MSE[c(2,4,6,8)]
linreg_cod = lr$test_MSE[c(1,3,5,7)]
linreg_cod_nt = lr$test_MSE[c(2,4,6,8)]

mses = data.frame( nn_cod, nn_cod_nt, linreg_cod, linreg_cod_nt)
mses = rbind(mses, c(NA, NA, struc$mean_MSE[1:2])) # hack for structure

# put all the MSEs in the right order to hack a stripchart
full2 = data.frame(t(full[,1:10]))
all = cbind( full2[,1:2], NA, NA, full2[,3:4], NA, NA, full2[,5:6], NA, NA, full2[,7:8], NA, NA, NA, NA, t(struc[1:10]))

label = c("A site", "-3:+2", "-5:+4", "-7:+5", "")
row.names(mses) = label
names(mses) = series_names

colors = c("mediumpurple1", "purple", "cornflowerblue", "darkslateblue")
struc.colors = c("firebrick1", "firebrick")
all.colors = c( rep(colors, 4), NA, NA, struc.colors)

pdf( out_fname, width = 2, height = 1.67, pointsize = 7, useDingbats = F, bg = "white" )
#cairo_pdf( out_fname, width = 2, height = 1.67, pointsize = 7 )
par( mex = 0.65 ) # sets margin stuff (stupidly)
par( mar = c(6,4.5,7,3) )
par( oma = c(0,1.5,1,0) ) # top and side margin for plot panel label
par( lwd = 0.75 )
par( xpd = NA )
centers = barplot( t(mses),
                   beside = T,
                   ylim = c( 0,1 ),
                   col = all.colors,
                   border = NA,
                   space = c(0.1,1),
                   axes = F,
                   legend.text = T,
                   args.legend = list( x = "topright", border = NA, bty ="n", cex = 0.6, inset = c(.5,-0.7)),
                   cex.names = 0.7,
                   ylab = "MSE",
                   xlab = NA
)
axis( 2, at = c(0,0.5,1), lwd = 0.75 )
mtext( "C", font = 2, line = -3, side = 3, outer = T, adj = 0 )

points( rep(as.numeric(centers), each = 10), as.matrix(all), 
        pch = 16,
        col = rgb(0.2,0.2,0.2,0.1),
        cex = .4)

# structure hacks
axis( 1, at = mean(centers[3:4,5]), labels = c("structure"), lwd=0, cex.axis=0.7)
legend( x = "topright", 
        legend = c("footprint","downstream"), 
        fill = struc.colors,
        inset = c(0, -0.7), 
        cex = 0.6,
        border = NA, bty="n")

abline( h = min(mses, na.rm=T), lty = 3, lwd = 0.5, xpd=F )

dev.off()

