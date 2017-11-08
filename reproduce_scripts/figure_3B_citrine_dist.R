# plot the distribution of citrine scores

args <- commandArgs(trailingOnly = TRUE)
random_citrine_scores <- args[1]
natural_scores <- args[2]
citrine_construct_scores <- args[3]
out_fname <- args[4]

sample = read.delim(random_citrine_scores,header=F)
natural = read.delim(natural_scores, header=T)
cit = read.delim( citrine_construct_scores, header=T, row.names = 1, 
		  comment.char="#")
names = c("MIN", "CHA2", "Y000", "Y333", "Y666", "Y999", "MAX")
nn.scores = cit[names,]

cols = c("magenta3", "purple2", "royalblue2", "green2", "yellow2", "orange2", "red2")
darkgrey50 <- rgb(t(col2rgb("darkgrey")[,1]), max = 255, alpha = 128)
purple30 <- rgb(t(col2rgb("purple")[,1]), max = 255, alpha = 76)
#cols = c("magenta", "red", "purple", "blue", "cyan", "green", "orange")

# not including CHA2
pdf(out_fname, width=2, height=1.167, pointsize=7, useDingbats = F, bg="white" )
#cairo_pdf(out_fname, width=2, height=1.167, pointsize=7 )
par( mex = 0.65 ) # sets margin stuff
par( mar =c(7,2.5,3,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
hist( sample$V1, 
      breaks = 30, 
      col = purple30,
      border = NA,
      xlim = c(150,400),
      freq = F,
      axes = F, 
      xlab = NA,
      ylab = NA,
      main = NA )
axis( 1 )
g = hist( natural$avg <- score * 238,
         plot = F,
         breaks = 30 )
plot(g, add = T, freq = F, col = darkgrey50, border = NA)
points(nn.scores[c(1,3:7)], rep(0.0025, length(nn.scores[c(1,3:7)])), col = cols[c(1,3:7)], pch = 20)
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
legend( "topright", legend = c("random eCitrine", "natural yeast genes"), fill = c(purple30,darkgrey50),
       border = NA, bty = "n", cex = 0.6, inset = c(0, -0.2) )
text(nn.scores[c(1,7)],c(max(h$density)/3, max(h$density)/3),labels = c("fastest","slowest"), col = c("magenta3", "red2"), cex = 0.7)
mtext( "B", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 
dev.off()


