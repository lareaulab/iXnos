# correlations of codon weights with tAI
# figure 2C
# data: Weinberg yeast
# codon properties source: ??

args <- commandArgs(trailingOnly = TRUE)
cod_scores_fname = args[1]
yeast_codon_props_fname = args[2]
out_fname = args[3]

prop = read.delim(yeast_codon_props_fname, header=T, stringsAsFactors = F)

wb = read.delim(cod_scores_fname, header=F, stringsAsFactors = F)

codons = sort( apply( expand.grid( c("A","C","G","T"), c("A","C","G","T"), c("A","C","G","T")), 1, paste, collapse = "" ))
pos = -7:5

row.names(wb) = codons
colnames(wb) = pos

label1 = c(-7,NA,-5,NA,-3,NA,"P",NA,1,NA,3,NA,5)
label2 = c(NA,-6,NA,-4,NA,"E",NA,"A",NA,2,NA,4,NA)

tai.cor = -cor( wb, prop$tAI, use = "pairwise.complete.obs", method = "spearman" )[,1]
tai.pval = apply(wb, 2, function(x){ cor.test(x, prop$tAI, use="pairwise.complete.obs", method="spearman")$p.value})

cols = rep("gray90", length(tai.pval))
cols[tai.pval < 0.05/length(tai.pval)] = "gray60"

pdf( out_fname, width=2, height=1.67, pointsize=7, useDingbats = F, bg = "white" )
#cairo_pdf( out_fname, width=2, height=1.67, pointsize=7 )
par( mex = 0.65 )
par( mar =c(6,5.5,5,3) )
par( oma = c(0,1.5,1,0) )
par( xpd = NA )
centers.tai = barplot( tai.cor,
                       xlab = "codon position",
                       ylab = "Spearman correlation",
                       col = cols,
                       border = cols,
                       space = 0.3,
                       axes = F,
                       axisnames = F )
axis( 2, seq( round(min(tai.cor)/2,1)*2, round(max(tai.cor/2),1)*2, by = 0.2 ), lwd = 0.75)
axis( 1, at = centers.tai, padj = -1,
      labels =  label1,
      tick=F, cex.axis = 0.7)
axis( 1, at = centers.tai, padj = -1,
      labels = label2,
      tick=F, cex.axis = 0.7)
mtext( "C", font = 2, line = -3, side = 3, outer = T, adj = 0 )
dev.off()
