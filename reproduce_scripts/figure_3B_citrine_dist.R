# plot the distribution of citrine scores

args <- commandArgs(trailingOnly = TRUE)
random_citrine_scores <- args[1]
citrine_construct_scores <- args[2]
out_fname <- args[3]

sample = read.delim(random_citrine_scores,header=F)
cit = read.delim( citrine_construct_scores, header=T, row.names = 1, 
		  comment.char="#")
names = c("CHA2", "MAX", "MIN", "Y000", "Y333", "Y666", "Y999")
nn.scores = cit[names,"nn.score"]

cols = c("magenta", "red", "purple", "blue", "cyan", "green", "orange")

# not including CHA2
pdf(out_fname, width=2, height=1.167, pointsize=7, useDingbats = F, bg="white" )
#cairo_pdf(out_fname, width=2, height=1.167, pointsize=7 )
par( mex = 0.65 ) # sets margin stuff
par( mar =c(6,2.5,3,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
hist( sample$V1, 
      breaks = 30, 
      col = "lightgrey",
      border = NA,
      xlim=c(150,400), 
      axes=F, 
      xlab = "NN score", 
      ylab = NA,
      main = NA )
axis( 1 )#, lwd = 0, lwd.ticks = 1 )
points(nn.scores[2:7], rep(750, length(nn.scores[2:7])), col = cols[2:7], pch = 20)
text(nn.scores[c(3,2)],c(10000,10000),labels = c("fastest","slowest"), col = c("purple", "red"))
mtext( "B", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 
dev.off()


