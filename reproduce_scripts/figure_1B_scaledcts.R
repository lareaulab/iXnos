# example of scaled counts from real gene
# figure 1B
# data: Weinberg yeast

args <- commandArgs(trailingOnly = TRUE)
id = args[1] # eg YOL086C for ADH1
g.range1 = as.numeric(args[2]) # eg 166
data_fname = args[3] # eg weinberg.te_set_bounds.size.27.31.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt
out_fname = args[4]

data = read.delim( data_fname, header = T)

# id = "YOL086C" # choose your gene
# g.range = 166:177 # choose your 12 codons
g.range = g.range1 : (g.range1 + 12 - 1)

g = data[data$gene == id,]

# make the coordinates work out with the missing 20 codons
prepend = data.frame(gene = rep(id, 20), 
                     cod_idx = 0:19, cod_seq = rep(NA, 20), 
                     raw_cts = rep(NA, 20), scaled_cts = rep(NA, 20))
g = rbind( prepend, g )

# calclulate the average coverage
xmax = max( g$cod_idx, na.rm = T )
ymean = sum( g$raw_cts, na.rm = T ) / (xmax - 20 + 1)

# the 12 codons we're plotting
codons = g$cod_seq[g.range]

# fit all the codon names in neatly as equally spaced nt
codons.array = unlist( strsplit( codons, split = ""))
label.pos = seq(1/6+0.1, 12, by = 1/3)

# put number tick marks only at multiples of 5
xticks = g.range
xticks[(xticks %% 5) != 0] = NA

# label and color the fast and slow codons
fast = which.min(g$raw_cts[g.range])
slow = which.max(g$raw_cts[g.range])
f = (fast - 1) * 3
s = (slow - 1) * 3
cols = rep( "black", 36)
cols[f+1:3] = "forestgreen"
cols[s+1:3] = "red"

# plot it
pdf( out_fname, width=2, height=1.67, pointsize=7, useDingbats = F, bg = "white" )
#cairo_pdf( out_fname, width=2, height=1.67, pointsize=7 )
par( mex = 0.65 ) # sets margin stuff (stupidly)
par( mar = c(6,4.5,7,4) )
par( oma = c(0,1.5,1,0) ) # top and side margin for plot panel label
par( lwd = 0.75 )
par( xpd = NA )
barplot(g$raw_cts[g.range],
        width = 0.8,
        space = 0.25, 
        border = NA, axes = F, 
        col="darkgray")

# average line
segments(x0 = 0.1, y0 = ymean, x1 = 12.1, y1 = ymean, lty=3)

# add all codons as nucleotides
text( label.pos, -700, labels = codons.array, cex = 0.4, font = 2, col = cols)

# label the fast and slow codons
arrows( x0 = fast-0.6, y0 = ymean, x1 = fast-0.6, y1 = g$raw_cts[g.range[fast]], col = "forestgreen", length = 0.02, code = 3)
arrows( x0 = slow-0.6, y0 = ymean, x1 = slow-0.6, y1 = g$raw_cts[g.range[slow]], col = "red", length = 0.02, code = 3)
text( fast-0.3, mean( c(ymean, g$raw_cts[g.range[fast]]) ), labels = c("fast"), col = "forestgreen", srt = 90, cex = 0.75)
text( slow-0.3, mean( c(ymean, g$raw_cts[g.range[slow]]) ), labels = c("slow"), col = "red", srt = 90, cex = 0.75)

# tick mark at all multiples of 5
axis(1, at = 0.6:11.6, labels = xticks, lwd=0, cex.axis = 0.7)

# raw and scaled count axes
axis(2, lwd=0.75, las=1, cex.axis = 0.5)
axis(4, lwd=0.75, las=1, at = c(0, ymean, 2*ymean, 3*ymean), labels = c(0,1,2,3))
title(ylab="raw count")
mtext( "scaled count", side = 4, line = 2)

# panel label
mtext( "B", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 
dev.off()

