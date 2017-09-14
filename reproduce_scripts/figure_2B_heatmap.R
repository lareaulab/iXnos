# heatmap of contributions of each codon at each position
# figure 2B
# data: Weinberg yeast

library("fields")

options("preferRaster" = F) # set just in case; later for scale bar will set to T

args <- commandArgs(trailingOnly = TRUE)
cod_scores_fname = args[1]
out_fname = args[2]

wb = read.delim(cod_scores_fname, header = F, stringsAsFactors = F, na.strings = "nan", colClasses = "numeric" )

codons = sort( apply( expand.grid( c("A","C","G","T"), c("A","C","G","T"), c("A","C","G","T")), 1, paste, collapse = "" ))

nonstop = intersect( intersect(which(codons != "TAG"), which(codons != "TAA")), which(codons != "TGA") )

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#AA0000"))

lastnt = rep(c("A","C","G","T"), 16)
firstnt = sort( apply( expand.grid( c("A","C","G","T"), c("A","C","G","T")), 1, paste, collapse = "" ))
firstnt = as.vector( rbind( firstnt, NA, NA, NA ))
firstnt[which(firstnt=="TA")+1] = "TA"
firstnt[which(firstnt=="TG")+1] = "TG"

cairo_pdf(out_fname, width=2, height=5, pointsize=7 )
par( mex = 0.65 )
par( mar =c(6,1.5,5,5) )
par( oma = c(0,1.5,1,0) )
image(-7:5,1:61,
      t(as.matrix(wb[rev(nonstop),])), 
      col = jet.colors(256), 
      asp=1, axes=F, 
      xlab="codon position", 
      ylab="")
# codon names (big letter for first two nt, small letter for third nt)
axis( 2, at = 1:61, labels = lastnt[rev(nonstop)], pos=-7, hadj=0.5, padj=0.5, las = 1, lwd = 0.75, tick = F, cex.axis = 0.6, family="mono" )
axis( 2, at = 1:61, labels = firstnt[rev(nonstop)], pos=-7.75, hadj=0.85, padj=0.7, las = 1, lwd = 0.75, tick = F, cex.axis = 1, family="mono" )

# group by first two nt, ticks / dividers sticking out on either side
axis( 2, at = c( 0.5, 4.5, 7.5, 11.5, seq( 13.5, 61.5, by = 4)), pos=-7.5, tcl = -2.5, lwd = 0.75, labels = F, tick = T )
axis( 2, at = c( 0.5, 4.5, 7.5, 11.5, seq( 13.5, 61.5, by = 4)), pos=5.5, tcl = 0.2, lwd = 0.75, labels = F, tick = T )

# label codon positions
label1 = c(-7,NA,-5,NA,-3,NA,"P",NA,1,NA,3,NA,5)
label2 = c(NA,-6,NA,-4,NA,"E",NA,"A",NA,2,NA,4,NA)
axis( 1, at = -7:5, labels = NA, tick = T, lwd = 0, lwd.ticks = 0.5)
axis( 1, at = seq(-7*1.5, 5*1.5, by=1.5), 
      labels = label1, 
      tick = F, lwd = 0, 
      cex.axis = 0.7 ) 
axis( 1, at = seq(-7*1.5, 5*1.5, by=1.5), 
      labels = label2, 
      tick = F, lwd = 0, 
      cex.axis = 0.7 ) 

# add top and bottom line since rect didn't dtrt
axis( 1, at = c(-7.5,5.5), labels = NA, tick = T, lwd.ticks=0, pos=0.5, lwd=0.75 )
axis( 1, at = c(-7.5,5.5), labels = NA, tick = T, lwd.ticks=0, pos=61.5, lwd=0.75 )

# label the legend (has to come before image.plot...)
mtext("faster", side=4, line = -2.5, adj = 0.265)
mtext("slower", side=4, line = -2.5, adj = 0.73)

mtext( "B", font = 2, line = -3, side = 3, outer = T, adj = 0 )

# diagonal lines
par(xpd=NA)
segments(  seq(-7*1.5, 5*1.5, by=1.5), -0.8,
           -7:5, -0.05,
           lwd = 0.5
)

options("preferRaster" = T) # rasterize the scale bar only

image.plot( zlim=c(min(wb,na.rm=T), max(wb,na.rm=T)),
            col = jet.colors(256),
            add = T, legend.only = T,
            legend.shrink = 0.5,
            legend.mar = 4,
            axis.args=list(las=0,lwd=0.75)
)
dev.off()
