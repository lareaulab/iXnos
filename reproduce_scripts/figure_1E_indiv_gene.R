# individual example gene
# figure 1E
# data: Weinberg yeast

args <- commandArgs(trailingOnly = TRUE)
id = args[1] # eg YOL086C for ADH1
data_fname = args[2] # eg weinberg.te_set_bounds.size.27.31.trunc.20.20.min_cts.200.min_cod.100.top.500.data_table.txt
y_te_hat_fname = args[3] # eg nonneg_str_epoch30_y_te_hat.csv
id_fname = args[4] #S_cer_id_symbol.txt
out_fname = args[5]

data = read.delim( data_fname, header = T)
y_hat = read.delim( y_te_hat_fname, header=F)
data$yhat = y_hat$V1

### gene name matching
namemap = read.delim( id_fname, header=F, row.names=1 )

g = data[data$gene == id,]

prepend = data.frame(gene = rep(id, 20), 
                     cod_idx = 0:19, cod_seq = rep(NA, 20), 
                     raw_cts = rep(NA, 20), scaled_cts = rep(NA, 20), 
                     yhat = rep(NA,20))

g = rbind( prepend, g )

symbol = namemap[id,1]
if(is.na(symbol)) { symbol = id }

ymax = max( g$scaled_cts, g$yhat, na.rm = T )
xmax = max( g$cod_idx, na.rm = T )
raw.ymax = max( g$raw_cts, na.rm = T )

ymean = sum( g$raw_cts, na.rm = T ) / (xmax - 20) # check these numbers...

pdf( out_fname, height=1.67, width=4, pointsize=7, useDingbats = F, bg = "white" )
par( mex = 0.65 ) 
par( mar =c(6,5.5,7,2) )
par( oma = c(0,1.5,1,0) ) # top and side margin for plot panel label
par( xpd = NA )
par( lwd = 0.75 )

bar = barplot(g$scaled_cts, 
              width=0.8, space=0.25, 
              border=NA, axes=F, 
              col="darkgray")
lines(bar, g$yhat, col="red",lwd=1)
title(xlab="codon position", ylab="scaled count")
axis(1, at = seq(0, xmax, by=100), lwd=0, lwd.ticks=0.75)
axis(2, at = seq(0, ymax, by=2), lwd = 0.75)#lwd=0, lwd.ticks=1, las=1)

mtext( paste(symbol), font = 3, line = 2, side = 3, outer = F, adj = 0 )
mtext( "E", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 

dev.off()
