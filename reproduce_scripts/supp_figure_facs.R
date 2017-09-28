args <- commandArgs(trailingOnly = TRUE)
facs_data_file = args[1]
citrine_construct_scores_file = args[2]
out_file = args[3]

facs.data = read.delim(facs_data_file, sep=",", header=T )
cit = read.delim(citrine_construct_scores_file, header=T, row.names = 1, comment.char="#")

medians = data.frame(tapply( facs.data$green / facs.data$red, list(facs.data$Isolate, facs.data$Strain), median ))
vars = data.frame(tapply( facs.data$green / facs.data$red, list(facs.data$Isolate, facs.data$Strain), var ))

medians = subset(medians, select = -c(NIY110, CHA1))
vars = subset(vars, select = -c(NIY110, CHA1))

bad = list( 
  c("CHA2", 4), # Nick called bad FACS data
  c("Y000", 8), # Nick called bad FACS data
  c("Y333", 1), # mutation
  c("Y333", 2), # Nick called bad FACS data
  c("Y333", 4), # mutation
  c("Y333", 7), # mutation
  c("Y333", 8), # mutation (status unknown)
  c("Y999", 3)  # copy number qPCR
)

lapply(bad, function(x) {
  medians[x[2], x[1]] <<- NA
  vars[x[2], x[1]] <<- NA
})

nn.scores = as.data.frame( sapply( names(medians), function(x) { rep( cit[x,"nn.score"], 8 )}) )

cols = rep( c("magenta", "red", "purple", "blue", "cyan", "green", "orange"), each = 8)

pdf(out_file, width=2, height=1.67, pointsize=7, useDingbats=F, bg="white" )
#cairo_pdf(out_file, width=2, height=1.67, pointsize=7 )
par( mex = 0.65 ) # sets margin stuff
par( mar =c(6,6.5,5,3) )
par( oma = c(0,0.5,1,0) )
plot(unlist(nn.scores), unlist(medians), 
     col = cols,
     axes = F,
     xlim = c(150,400),
     ylim = c(0, max(medians,na.rm=T)),
     cex = 0.6,
     pch = 3,
     #     pch = 20,
     xlab = "NN score",
     ylab = "eCitrine / mCherry\nfluorescence ratio"
)
axis( 1 )
axis( 2, at = seq(0, max(medians,na.rm=T), by=0.5))
mtext( "A", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 
dev.off()

