
codons = sort( apply( expand.grid( c("A","C","G","T"), c("A","C","G","T"), c("A","C","G","T")), 1, paste, collapse = "" ))
pos = -7:5

cl2 = read.delim("lareau_s28_codon_scores.tsv", header=F, stringsAsFactors = F, colClasses = "numeric", row.names = codons, col.names = pos, na.strings="nan")
cl1 = read.delim("green_s28_codon_scores.tsv", header=F, stringsAsFactors = F, colClasses = "numeric", row.names = codons, col.names = pos, na.strings="nan")

tested = c("ATA", "TCC", "CCA", "CGT", "GAC", "GGG")

scores = data.frame( cl1.green = round( cl1[tested,]$X.5, 3 ), 
                     cl2.lareau = round( cl2[tested,]$X.5, 3 )
)
row.names(scores) = tested


cairo_pdf("supp_figure_cl1v2.pdf", width=3, height=3, pointsize=7 )
par( mex = 0.65 ) # sets margin stuff
par( mar =c(6,6.5,5,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot( cl2$X.5, 
      cl1$X.5, 
      xlab = "5' codon weights, circligase II",
      ylab = "5' codon weights, circligase I",
      pch = 20,
      bty = "n",
      col = "darkgrey"
      )

points( cl2[tested,]$X.5,
        cl1[tested,]$X.5,
        col="red",
        pch=20 )

text( cl2[tested,]$X.5,
      cl1[tested,]$X.5,
      labels = tested,
      pos = 1 )
dev.off()
