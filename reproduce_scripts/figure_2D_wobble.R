# wobble base pairing in P site
# figure 2D
# data: weinberg yeast
# codon: anticodon source: compiled from Johansson et al, MCB, 2008

args <- commandArgs(trailingOnly = TRUE)
cod_scores_fname = args[1] # weinberg_rep3_L1W200_ff_n7p5_n21p17_rep0_feat_imps.tsv or equiv.
codon_fname = args[2] # yeast_codon_anticodon.csv
out_fname = args[3]

bp = read.delim(codon_fname, sep=",", header=T)
# to make the codon_anticodon file, we chose the direct-paired anticodon if there were two anticodon options, or the standard wobble pair in the case of CCC.

codons = sort( apply( expand.grid( c("A","C","G","U"), c("A","C","G","U"), c("A","C","G","U")), 1, paste, collapse = "" ))
pos = -7:5

wb = read.delim( cod_scores_fname,
                 header = F, stringsAsFactors = F,
                 na.strings = "nan", colClasses = "numeric",
                 row.names = codons, col.names = pos
                 )

bp$P.site.score = wb$X.1[match(bp$codon, row.names(wb))]

s.cols = c(
 "I:C" = "green",
 "I:A" = "red",
 "G:U" = "purple",
 "I:U" = "blue",
 "G:C" = "darkgrey",
 "U:A" = "darkgrey"
)

bp$pair = bp$simplified
bp$pair[bp$simplified == "C:G"] = "G:C"
bp$pair[bp$simplified == "U:G"] = "G:U"

mylevels = levels(bp$pair)
mylevels[mylevels == "C:G"] = "G:C"
mylevels[mylevels == "U:G"] = "G:U"
levels(bp$pair) = mylevels
## put the stripchart in order of average score
myorder = order(tapply( bp$P.site.score, bp$pair, mean, simplify = T))

## figure out if any columns are significant
len = dim(bp)[1]
gc = which(bp$pair == "G:C")
ua = which(bp$pair == "U:A")
ic = which(bp$pair == "I:C")
iu = which(bp$pair == "I:U")
gu = which(bp$pair == "G:U")
ia = which(bp$pair == "I:A")

P.pvals = list()
P.pvals["I:C"] = wilcox.test(bp$P.site.score[ic], bp$P.site.score[setdiff(1:len, ic)], alternative="t", paired=F)$p.value
P.pvals["I:U"] = wilcox.test(bp$P.site.score[iu], bp$P.site.score[setdiff(1:len, iu)], alternative="t", paired=F)$p.value
P.pvals["I:A"] = wilcox.test(bp$P.site.score[ia], bp$P.site.score[setdiff(1:len, ia)], alternative="t", paired=F)$p.value
P.pvals["G:U"] = wilcox.test(bp$P.site.score[gu], bp$P.site.score[setdiff(1:len, gu)], alternative="t", paired=F)$p.value
P.pvals["G:C"] = wilcox.test(bp$P.site.score[gc], bp$P.site.score[setdiff(1:len, gc)], alternative="t", paired=F)$p.value
P.pvals["U:A"] = wilcox.test(bp$P.site.score[ua], bp$P.site.score[setdiff(1:len, ua)], alternative="t", paired=F)$p.value
P.pvals = unlist(P.pvals[levels(bp$pair)[myorder]])

n = length(levels(bp$pair))

cairo_pdf( out_fname, width=2, height=1.67, pointsize=7)
par( mex = 0.65 ) 
par( mar =c(6,5.5,5,3) )
par( oma = c(0,1.5,1,0) )
par( lwd = 0.75 )
stripchart( P.site.score ~ factor(pair, levels = levels(pair)[myorder]), 
            data = bp,
            method = "jitter", 
            jitter = 0.15, 
            vertical = T, 
            pch = 20, cex = 0.5,
            frame.plot = F,
            axes = F,
            ylab = "P site weight",
            col = s.cols[levels(bp$pair)[myorder]])
axis( 2, at = seq(0, 1.2, by = 0.6), lwd = 0.75)
axis( 1, at = 1:n, labels = levels(bp$pair)[myorder], lwd = 0, cex.axis = 0.7)
axis( 1, pos = -0.3, lwd = 0,
      at = which( P.pvals * n < 0.05 ), labels = rep("**",length(which(P.pvals * n < 0.05))))

mtext( "D", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 
dev.off()

