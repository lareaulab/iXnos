# comparison of three methods
# supp figure, panel B
# Liu data from running riboshape on the Weinberg dataset locally
# O'Connor data from running RUST on 28mers from the Weinberg dataset via their Galaxy server (offset = 15)

args <- commandArgs(trailingOnly = TRUE)
tunney_fname = args[1] # corrs_by_gene_density.txt
oconnor_fname = args[2] # oconnor_rust_weinberg.csv
liu_fname = args[3] # ASK: corrs_by_gene.txt
out_fname = args[6] 

# our data
tunney = read.delim( tunney_fname, header=T, comment.char = "#", row.names=1 )

tunney.names = row.names(tunney)

num = length(tunney.names)
training = which( as.logical(tunney$tr_gene) )
nottraining = which( !as.logical(tunney$tr_gene) )

oconnor = read.delim( oconnor_fname, header=T, sep=",", row.names=1 )

liu = read.delim( liu_fname, header=F, row.names=1)

#### combine them all -- genes that are found in all
in.all = intersect( row.names(oconnor), intersect( tunney.names[nottraining], row.names(liu) ))

ord = order( tunney[in.all,"fp_density"], decreasing=T )
in.all = in.all[ord] # sort gene names by density

num.all = length(in.all)

oconnor.x = match( in.all, row.names(oconnor) )
liu.x = match( in.all, row.names(liu) )
tunney.x = match( in.all, row.names(tunney) )

tunney.loess = predict( loess( tunney$pearson_r[tunney.x] ~ c(1:num.all) ),
                        1:num.all, se=T )
oconnor.loess = predict( loess( oconnor$Pearson.s.coefficient[oconnor.x] ~ c(1:num.all) ),
                         1:num.all, se=T )
liu.loess =  predict( loess( liu[liu.x,8] ~ c(1:num.all) ),
                      1:num.all, se=T )

all.data = data.frame( t = tunney$pearson_r[tunney.x],
                       o = oconnor$Pearson.s.coefficient[oconnor.x],
                       l = liu[liu.x, 8]
)

loess.fits = data.frame( t = tunney.loess$fit,
                         o = oconnor.loess$fit,
                         l = liu.loess$fit
)

loess.ses = data.frame( t = tunney.loess$se.fit,
                        o = oconnor.loess$se.fit,
                        l = liu.loess$se.fit
)

cols = c("red", "blue", "green")
titles = c("Tunney", "O'Connor", "Liu")

make.plot = function( x, y, z, w) {
  plot( 1:num.all,
        x,
        pch = 21, 
        cex = 0.5,
        col = NA,
        bg = "#00000040",
        bty = "n",
        axes = F,
        ylim = c(0, 1),
        xlab = NA,
        ylab = NA
  )
  lines( 1:num.all, y, col=z, lwd = 1.5 )
  axis( 2, at = c(0,0.5,1), labels = c(0, NA, 1), las = 1)
  mtext( w, side = 3, line = -1.5, adj = 1, col=z )
  
}

# 99% confidence intervals
se.polygon = function( x, y ){ 
  y.polygon <- c( (x + 2.576 * y), rev(x - 2.576 * y) )
  x.polygon <- c( 1:num.all, num.all:1 )
  polygon(x.polygon, y.polygon, col="#00009933", border=NA)
}

pdf( out_fname, width=2, height=3.5, pointsize=7, useDingbats=F )
par( mfrow = c(4,1) )
par( mex = 0.65 ) # sets margin stuff
par( lwd = 0.75 )
par( xpd = NA )
par(oma = c(5, 1.5, 6, 0),
    mar = c(1, 5.5, 0, 3) )
par( cex = 1 )

Map( make.plot, all.data, loess.fits, cols, titles )

## compare loess lines of all
plot( NA,
      xlim = c(1, num.all),
      ylim = c(0,1),
      bty = "n",
      axes = F,
      xlab = NA,
      ylab = NA
)

# plot SE margins
Map( se.polygon, loess.fits, loess.ses )
# plot all loess lines
Map( function(x,y){ lines( 1:num.all, x, col=y) }, loess.fits, cols )
axis( 1 )
axis( 2, at = c(0,0.5,1), labels = c(0, NA, 1), las = 1)
title( xlab = "genes sorted by footprint density" )
mtext( "Pearson correlation per gene", side = 2, outer = T, line = -2.5)
mtext( "B", font = 2, line = 2, side = 3, outer = T, adj = 0 ) 
dev.off()

# RUST28
# > sapply(all.data, mean)
# t         o         l       l.s 
# 0.5590593 0.4802878 0.4108693 0.4750529 


