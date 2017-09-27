# plot showing A site weights in yeast and human data generated with same circligase
# figure 2G
# data: Lareau / Graham yeast, Iwasaki human

args <- commandArgs(trailingOnly = TRUE)
lareau_cod_scores_fname = args[1]
iwasaki_cod_scores_fname = args[2]
out_fname = args[3]

codons = sort( apply( expand.grid( c("A","C","G","T"), c("A","C","G","T"), c("A","C","G","T")), 1, paste, collapse = "" ))
pos = -7:5

iwasaki = read.delim(iwasaki_cod_scores_fname, header=F, stringsAsFactors = F, colClasses = "numeric", row.names = codons, col.names = pos, na.strings="nan")
lareau = read.delim(lareau_cod_scores_fname, header=F, stringsAsFactors = F, colClasses = "numeric", row.names = codons, col.names = pos, na.strings="nan")

a.cor = round( cor(iwasaki$X0, lareau$X0, method="spearman", use="complete.obs"), 2)
end.cor = round( cor(iwasaki$X.5, lareau$X.5, method="spearman", use="complete.obs"), 2)

xmin = min( iwasaki$X0, na.rm = T )
xmin = min( xmin, round(xmin) )
xmax = max( iwasaki$X0, na.rm = T )
xmax = max( xmin, round(xmax) )
xat = round(xmin):round(xmax)
xlab = xat
xlab[ which(xat %% 2 != 0) ] = NA

ymin = min( lareau$X0, na.rm = T )
ymin = min( ymin, round(ymin) )
ymax = max( lareau$X0, na.rm = T )
ymax = max( ymin, round(ymax) )

#pdf( out_fname, width=2, height=1.67, pointsize=7, useDingbats = F, bg = "white" )
cairo_pdf( out_fname, width=2, height=1.67, pointsize=7 )
par( mex = 0.65 )
par( mar = c(6,8.5,5,6) )
par( oma = c(0,1.5,1,0) )
par( lwd = 0.75 )
par( xpd = NA )
plot( iwasaki$X0, lareau$X0, 
      axes = F,
      cex = 0.5,
      xlim = c( xmin, xmax ),
      ylim = c( ymin, ymax ),
      pch = 20, col="grey50", 
      xlab = "A site, human", 
      ylab = "A site, yeast")
axis( 1, lwd = 0.75, at = xat, labels = xlab )
axis( 2, lwd = 0.75, at = round(ymin):round(ymax) )
mtext( "G", font = 2, line = -3, side = 3, outer = T, adj = 0 )
mtext( bquote(rho==.(a.cor)), side = 1, line = -2, adj = 1 )
dev.off()
