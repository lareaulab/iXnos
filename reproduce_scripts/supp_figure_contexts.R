## supp figure to identify contexts that do better than their A site
## data: Weinberg A site model and -3:2 model with no A site

args <- commandArgs(trailingOnly = TRUE)
y_fname = args[1] # /mnt/lareaulab/rtunney/iXnos/expts/weinberg/process/te_set_bounds.size.27.31.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt
# asite_fname = args[2] # /mnt/lareaulab/rtunney/iXnos/expts/weinberg/lasagne_nn/full_cod_p0_nt_p0p2_rep0/epoch60/y_te_hat.txt
# noa_fname = args[3] # /mnt/lareaulab/rtunney/iXnos/expts/weinberg/lasagne_nn/noAsite_cod_n3p2_nt_n9p8/epoch60/y_te_hat.txt
# out_list_fname = args[4] # high_asite.csv - write the high A site error output
contexts_fname = args[5] # /mnt/lareaulab/rtunney/iXnos/results/paper_data/high_asite_contexts.txt - pull in the contexts from those sites (separate perl script...)
out_fname = args[6] # supp_figure_contexts.pdf

# observed data
y = read.delim( y_fname, header=T )
# # A site only
# asite = read.delim( asite_fname, header=F )
# # no a site, n3p2 context
# noa = read.delim( noa_fname, header=F )
# 
# asite.error = (y$scaled_cts - asite[,1])^2
# noa.error = (y$scaled_cts - noa[,1])^2
# 
# # more error is due to A site, for each codon (high = other stuff is a better predictor)
# high.asite = asite.error > noa.error
#   
# write.table(y[high.asite, 1:3], out_list_fname, sep=",", quote=F, row.names = F)
# 
# # then run the program that pulls out their contexts
# # ./get_contexts.pl > high_asite_contexts.txt

bg_codon_freqs = as.data.frame( table(y$cod_seq), row.names = 1 )

contexts = read.delim( contexts_fname, sep = " ", header=F)
names(contexts) = c("cod-3", "codE", "codP", "codA", "cod1", "cod2")

codoncounts = as.data.frame(sapply( names(contexts), function(x){ table(contexts[,x]) }))

n.all = sum(bg_codon_freqs$Freq)
n.high_asite = sum(codoncounts$codA)

all = expand.grid( row.names(bg_codon_freqs), names(codoncounts) )

prop.test.all = apply( all, 1, function(x) {
  prop.test( x = c( codoncounts[x[1],x[2]], bg_codon_freqs[x[1],1] ),
             n = c( n.high_asite, n.all ),
             alternative = "two.sided"
  )$p.value
})

prop.test.all =  matrix( prop.test.all, nrow = 61, ncol = 6) 
row.names(prop.test.all) = row.names(bg_codon_freqs)
colnames(prop.test.all) = names(codoncounts)

prop.test.2 = matrix( p.adjust(prop.test.all, method = "fdr"), nrow = 61, ncol = 6 )
row.names(prop.test.2) = row.names(bg_codon_freqs)
colnames(prop.test.2) = names(codoncounts)

prop.test.2[,4] = NA # A site
prop.test.2[ prop.test.2 > 0.05 ] = NA

codonrange1 = -3
codonrange2 = 2
label1 = c(-3, NA,  "P", NA,  "+1",  NA)
label2 = c(NA, "E", NA,  NA, NA, "+2") # don't include A site

codons = sort( apply( expand.grid( c("A","C","G","T"), c("A","C","G","T"), c("A","C","G","T")), 1, paste, collapse = "" ))
nonstop = intersect( intersect(which(codons != "TAG"), which(codons != "TAA")), which(codons != "TGA") )

lastnt = rep(c("A","C","G","T"), 16)
firstnt = sort( apply( expand.grid( c("A","C","G","T"), c("A","C","G","T")), 1, paste, collapse = "" ))
firstnt = as.vector( rbind( NA, NA, firstnt, NA ))

firstnt2 = rep(NA, 64)
firstnt2[which(firstnt=="TA")-1] = "TA"

### plot the heatmap
pdf( out_fname, width=1.67, height=5, pointsize=7, useDingbats = F, bg = "white" )
par( mex = 0.65 )
par( mar =c(4.5,5,5,6) )
par( oma = c(1.5,1,0,0) )
# par( mar =c(6,1.5,5,5) )
# par( oma = c(0,1.5,1,0) )
image( codonrange1:codonrange2, 1:61,
       t(prop.test.2[61:1,]), 
       asp = 1, axes = F, 
       zlim = c(0, 0.05),
       xlab = NA,
       ylab = NA ,
       useRaster = F
       )
# codon names (big letter for first two nt, small letter for third nt)
axis( 2, at = 1:61, labels = lastnt[rev(nonstop)], pos = codonrange1 + 0.2, hadj=0.5, padj=0.5, las = 0, lwd = 0.75, tick = F, cex.axis = 0.6, family="mono" )
axis( 2, at = 1:61, labels = firstnt[rev(nonstop)], pos = codonrange1 - 1, hadj=0.85, padj=0.7, las = 0, lwd = 0.75, tick = F, cex.axis = 0.9, family="mono" )
axis( 2, at = 1:61, labels = firstnt2[rev(nonstop)], pos = codonrange1 - 1, hadj=0.85, padj=0.7, las = 0, lwd = 0.75, tick = F, cex.axis = 0.9, family="mono" )

# group by first two nt, ticks / dividers sticking out on either side
axis( 2, at = c( 0.5, 4.5, 7.5, 11.5, seq( 13.5, 61.5, by = 4)), pos = codonrange1 - 0.5, tcl = -2.5, lwd = 0.75, labels = F, tick = T )
axis( 2, at = c( 0.5, 4.5, 7.5, 11.5, seq( 13.5, 61.5, by = 4)), pos = codonrange2 + 0.5, tcl = 0.2, lwd = 0.75, labels = F, tick = T )

# label codon positions
axis( 1, at = codonrange1:-1, labels = NA, tick = T, lwd = 0, lwd.ticks = 0.5)
axis( 1, at = 1:codonrange2, labels = NA, tick = T, lwd = 0, lwd.ticks = 0.5)
axis( 1, at = seq(codonrange1*1.5, codonrange2*1.5, by=1.5), 
      labels = label1, 
      tick = F, lwd = 0,
      line = 1.1,
      las = 2,
      cex.axis = 0.7 ) 
axis( 1, at = seq(codonrange1*1.5, codonrange2*1.5, by=1.5), 
      labels = label2, 
      tick = F, lwd = 0, 
      line = 1.1,
      las = 2,
      cex.axis = 0.7 ) 
# add top and bottom line since rect didn't dtrt
axis( 1, at = c( codonrange1 - 0.5, codonrange2 + 0.5 ), labels = NA, tick = T, lwd.ticks=0, pos=0.5, lwd=0.75 )
axis( 1, at = c( codonrange1 - 0.5, codonrange2 + 0.5 ), labels = NA, tick = T, lwd.ticks=0, pos=61.5, lwd=0.75 )

# diagonal lines
par(xpd=NA)
# don't do A site
segments(  seq(codonrange1*1.5, -1*1.5, by=1.5), -0.8,
           codonrange1:-1, -0.05,
           lwd = 0.5 )
segments(  seq(1*1.5, codonrange2*1.5, by=1.5), -0.8,
           1:codonrange2, -0.05,
           lwd = 0.5 )

options("preferRaster" = T) # rasterize the scale bar only

image.plot( zlim=c(0, 0.05),
            col = heat.colors(256),
            add = T, legend.only = T,
            legend.shrink = 0.5,
            legend.mar = 5,
            axis.args=list(las=0,lwd=0.75)
)
dev.off()

