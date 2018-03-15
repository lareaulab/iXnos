## supp figure to identify A sites that dominate the predictions over their surroundings
## data: Weinberg -3:2 model and -3:2 model with no A site

args <- commandArgs(trailingOnly = TRUE)
y_fname = args[1] # /mnt/lareaulab/rtunney/iXnos/expts/weinberg/process/te_set_bounds.size.27.31.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt
n3p2_fname = args[2] # /mnt/lareaulab/rtunney/iXnos/expts/weinberg/lasagne_nn/full_cod_n3p2_nt_n9p8_rep0/epoch60/y_te_hat.txt
n3p2noA_fname = args[3] # "/mnt/lareaulab/rtunney/iXnos/expts/weinberg/lasagne_nn/noAsite_cod_n3p2_nt_n9p8/epoch60/y_te_hat.txt"
out_fname = args[4] # supp_figure_a_site.pdf

# real scaled counts
y = read.delim( y_fname, header=T )
# no a site, n3p2 context
noa = read.delim( n3p2noA_fname, header=F)
# -3 to +2 including A site
n3p2 = read.delim( n3p2_fname, header=F)

all.cors = function( x, y, g ) { 
  cors <- sapply( levels(g), function(i) { cor(x[g == i], y[g == i]) } )
  names( cors ) <- levels(g)
  cors
}
cor.noa.percodon = all.cors( noa$V1, y$scaled_cts, y$cod_seq )
cor.n3p2.percodon = all.cors( n3p2$V1, y$scaled_cts, y$cod_seq )

counts = sapply( levels(y$cod_seq), function(i) { length(noa$V1[ y$cod_seq == i]) } )

fisher = function( r1, r2, n ) { ( atanh(r1) - atanh(r2) ) / (2/(n-3))^0.5 }

fisher.out = fisher( cor.noa.percodon, cor.n3p2.percodon, counts )
fisher.p = 2 * (1 - pnorm( abs(fisher.out) ))

sig = which( p.adjust(fisher.p, method="fdr") < 0.05 )


pdf( out_fname, width=2, height=1.67, pointsize=7, useDingbats=F)
par( mex = 0.65 ) # sets margin stuff
par( mar = c(6,6.5,4,3) )
par( oma = c(0,1.5,1,0))
par( lwd = 0.75 )
plot( cor.n3p2.percodon, 
      cor.noa.percodon, 
      xlab = "correlation, -3 to +2 region",
      ylab = "correlation, -3 to +2 region\nwithout A site",
      bty = "n",
      xaxs = "i",
      asp = 1,
      xlim = c(0,0.8),
      ylim = c(0,0.8),
      pch = 19,
      cex = 0.5
)
points( cor.n3p2.percodon[sig], cor.noa.percodon[sig], pch = 19, cex = 0.5, col = "red" )
text( cor.n3p2.percodon[sig], cor.noa.percodon[sig], labels = names(sig), pos = c(1,4), cex = 0.5, col = "red" )
abline( 0, 1, col = "grey", lty = 2 )
dev.off()
