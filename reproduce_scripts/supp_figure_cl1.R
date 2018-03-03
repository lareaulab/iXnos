# circligase I test
# supp figure
# data: Lareau yeast, Green yeast

args <- commandArgs(trailingOnly = TRUE)
lareau_scores_fname = args[1] # 28mers only!
green_scores_fname = args[2] # 28mers only!
qpcr_fname = args[3]
out_fname = args[4]

codonrange1 = -5
codonrange2 = 4

oligos = list( "1" = "ATA", "2" = "TCC", "3" = "CCA", "7" = "CGT", "8" = "GAC", "9" = "GGG")

codons = sort( apply( expand.grid( c("A","C","G","T"), c("A","C","G","T"), c("A","C","G","T")), 1, paste, collapse = "" ))
pos = codonrange1:codonrange2

cl2 = read.delim( lareau_scores_fname, header=F, stringsAsFactors = F, colClasses = "numeric", row.names = codons, col.names = pos, na.strings="nan")
cl1 = read.delim( green_scores_fname, header=F, stringsAsFactors = F, colClasses = "numeric", row.names = codons, col.names = pos, na.strings="nan")

tested = c("ATA", "TCC", "CCA", "CGT", "GAC", "GGG")

scores = data.frame( cl1.green = round( cl1[tested,]$X.5, 3 ), 
                     cl2.lareau = round( cl2[tested,]$X.5, 3 )
)
row.names(scores) = tested

qdat = read.delim( qpcr_fname, header=T, sep=",")
qdat = droplevels( qdat[ qdat$Oligo != "804", ] )

circprim = "NM827_NM828"
contprim = "NM828_804_LC_1"

se = function(x) { sd(x)/sqrt(length(x)) }

## CALCULATE EFFICIENCY AVERAGES
eff = tapply(qdat$eff, qdat$PRIMERS, mean)
eff.se = tapply(qdat$eff, qdat$PRIMERS, se)

cl1.cont = qdat[ qdat$PRIMERS == contprim & qdat$CL == 1, ]
cl1.circ = qdat[ qdat$PRIMERS == circprim & qdat$CL == 1, ]
cl2.cont = qdat[ qdat$PRIMERS == contprim & qdat$CL == 2, ]
cl2.circ = qdat[ qdat$PRIMERS == circprim & qdat$CL == 2, ]

## amplification amounts
cl1.circ.amp = eff[circprim]^cl1.circ$Cy0
cl2.circ.amp = eff[circprim]^cl2.circ$Cy0
cl1.cont.amp = eff[contprim]^cl1.cont$Cy0
cl2.cont.amp = eff[contprim]^cl2.cont$Cy0

cl1.circ.mean = tapply( cl1.circ.amp, cl1.circ$Oligo, mean)
cl2.circ.mean = tapply( cl2.circ.amp, cl2.circ$Oligo, mean)
cl1.cont.mean = tapply( cl1.cont.amp, cl1.cont$Oligo, mean)
cl2.cont.mean = tapply( cl2.cont.amp, cl2.cont$Oligo, mean)

names(cl1.circ.mean) = unlist(oligos[names(cl1.circ.mean)]) # so we can sanity-check that scores and results are paired properly
names(cl2.circ.mean) = unlist(oligos[names(cl2.circ.mean)])
names(cl1.cont.mean) = unlist(oligos[names(cl1.cont.mean)])
names(cl2.cont.mean) = unlist(oligos[names(cl2.cont.mean)])

# ratios
cl1.ratio = cl1.cont.mean/cl1.circ.mean
cl2.ratio = cl2.cont.mean/cl2.circ.mean
ratio.max = max( cl1.ratio, cl2.ratio )

cl1.ratio = cl1.ratio/ratio.max
cl2.ratio = cl2.ratio/ratio.max

# error propagation
cl1.circ.se = tapply( cl1.circ.amp, cl1.circ$Oligo, se)
cl2.circ.se = tapply( cl2.circ.amp, cl2.circ$Oligo, se)
cl1.cont.se = tapply( cl1.cont.amp, cl1.cont$Oligo, se)
cl2.cont.se = tapply( cl2.cont.amp, cl2.cont$Oligo, se)

cl1.se.percent = sqrt( (cl1.circ.se / cl1.circ.mean)^2 + (cl1.cont.se / cl1.cont.mean)^2 )
cl1.se.abs = cl1.ratio * cl1.se.percent

cl2.se.percent = sqrt( (cl2.circ.se / cl2.circ.mean)^2 + (cl2.cont.se / cl2.cont.mean)^2 )
cl2.se.abs = cl2.ratio * cl2.se.percent

ymax = max( cl1.ratio + cl1.se.abs, cl2.ratio + cl2.se.abs )

ord = order(scores$cl1.green)
down = ord[c(2,5)]
up = ord[c(1,3,4,6)]

pdf( out_fname, width=3, height=3, pointsize=7, useDingbats=F, bg="white" )
#pdf( out_fname, width=2, height=1.67, pointsize=7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) # sets margin stuff
par( oma = c(0,1.5,1,0) )
par( mar = c(6,6,5,6) )
par( lwd = 0.75 )
plot( scores$cl1.green, cl1.ratio,
      xlab = "bias score", ylab = "relative ligation", 
      xlim = c(-1, 1), 
      ylim = c(0,ymax),
      pch = NA, bty = "n"
)
abline( lm(cl1.ratio ~ scores$cl1.green), col = "darkgray", lty = 3 )
arrows( scores$cl1.green, cl1.ratio - cl1.se.abs,
        scores$cl1.green, cl1.ratio + cl1.se.abs,
        angle = 90, code = 3, length = 0.025, lwd = 0.75, col = "gray50" )
points( scores$cl1.green, cl1.ratio, pch=20, col = "blue" )
par( xpd = NA )
text( scores$cl1.green[down], cl1.ratio[down] - cl1.se.abs[down] - 0.075, labels = row.names(scores)[down], cex = 0.5, col = "gray50" )
text( scores$cl1.green[up], cl1.ratio[up] + cl1.se.abs[up] + 0.075, labels = row.names(scores)[up], cex = 0.5, col = "gray50" )
mtext( "B", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 
dev.off()
