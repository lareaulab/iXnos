# leave-one-out analysis of codon position contributions
# supp figure A
# data: Iwasaki human

args <- commandArgs(trailingOnly = TRUE)
leaveout_fname = args[1] # results/iwasaki/leaveout_series/leaveout_corrs.txt
full_fname = args[2] # results/iwasaki/feat_neighborhood_series/feat_neighborhood_corrs.txt
out_fname = args[3] # supp_figure_iwasakicorrs.pdf

full.name = "full_cod_n7p5_nt_n21p17"

dt = read.delim(leaveout_fname, header = T, comment = "#")
full = read.delim(full_fname, header = T, comment = "#")

leaveout.mean = apply(dt[1:13,2:11], 1, mean)

full.row = which( full$rep_series == full.name)
full.mean = apply(full[full.row,2:11], 1, mean)

data = data.frame( cbind( t(full[full.row,2:11]), NA, t(dt[1:13,2:11]) ))

label1 = c("all", NA, -7, NA, -5, NA, -3, NA,  "P", NA,  1,  NA, 3,  NA, 5)
label2 = c(NA,     NA, NA, -6, NA, -4, NA, "E", NA,  "A", NA, 2,  NA, 4,  NA)

pdf( out_fname, width=2, height=1.67, pointsize = 7, useDingbats = F, bg = "white" )
#cairo_pdf( out_fname, width=2, height=1.67, pointsize = 7 )
par( mex = 0.65 )
par( mar = c(6,5.5,4,3) )
par( oma = c(0,1.5,1,0) )
par( lwd = 0.75 )
par( xpd = NA )

plot( NA, 
      xlim = c( 1, 15 ),
      ylim = c( 0, 0.15 ),
      axes = F,
      xlab = "codon position",
      ylab = expression(paste(Delta, " correlation")))

rect( c(2.6:14.6), 0, c(3.4:15.4), full.mean - leaveout.mean,
      col = "mediumpurple3", border = NA )

stripchart( full.mean - data, 
            vertical = T, 
            pch = 16,
            col = rgb(0.2,0.2,0.2,.3),
            cex = .4,
            add = T
            )

axis( 2, at = seq(0,0.15,by=0.05), labels = c("0.0",NA,"0.1",NA), lwd = 0.75 )
axis( 1, at = 1:15, padj = -1, labels = label1, tick = F, cex.axis = 0.7)
axis( 1, at = 1:15, padj = -1, labels = label2, tick = F, cex.axis = 0.7)

mtext( "A", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 
dev.off()
