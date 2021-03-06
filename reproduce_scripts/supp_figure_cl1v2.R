# circligase I vs II codon scores
# supp figure 
# data: Lareau yeast, Green yeast

args <- commandArgs(trailingOnly = TRUE)
lareau_scores_fname = args[1] # 28mers only!
green_scores_fname = args[2] # 28mers only!
out_fname = args[3]

codonrange1 = -5
codonrange2 = 4

codons = sort( apply( expand.grid( c("A","C","G","T"), c("A","C","G","T"), c("A","C","G","T")), 1, paste, collapse = "" ))
pos = codonrange1:codonrange2

cl2 = read.delim( lareau_scores_fname, header=F, stringsAsFactors = F, colClasses = "numeric", row.names = codons, col.names = pos, na.strings="nan" )
cl1 = read.delim( green_scores_fname, header=F, stringsAsFactors = F, colClasses = "numeric", row.names = codons, col.names = pos, na.strings="nan" )

tested = c("ATA", "TCC", "CCA", "CGT", "GAC", "GGG")

scores = data.frame( cl1.green = round( cl1[tested,]$X.5, 3 ), 
                     cl2.lareau = round( cl2[tested,]$X.5, 3 )
)
row.names(scores) = tested

#pdf( out_fname, width=2, height=1.67, pointsize=7, useDingbats = F, bg = "white" )
pdf( out_fname, width=3, height=3, pointsize=7, useDingbats=F, bg="white" )
par( mex = 0.65 ) # sets margin stuff
par( lwd = 0.75 )
par( oma = c(0,1.5,1,0) )
par( mar = c(6,6,5,6) )
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

text( cl2[tested[2:6],]$X.5,
      cl1[tested[2:6],]$X.5,
      labels = tested[2:6],
      pos = 1 )
text( cl2[tested[1],]$X.5,
      cl1[tested[1],]$X.5,
      labels = tested[1],
      pos = 3 )


mtext( "A", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 

dev.off()
