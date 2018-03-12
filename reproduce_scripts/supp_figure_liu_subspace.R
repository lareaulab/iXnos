# analysis of Liu subspace performance
# Liu data from running riboshape on the Weinberg dataset locally

args <- commandArgs(trailingOnly = TRUE)
tunney_fname = args[1] # corrs_by_gene_density.txt
liu_s_fname = args[4] # ASK: subspace_corrs.txt
liu_gene_fname = args[5] # ASK: gene_names.txt
out_fname = args[6] 

# our data
tunney = read.delim( tunney_fname, header=T, comment.char = "#", row.names=1 )

tunney.names = row.names(tunney)

num = length(tunney.names)
training = which( as.logical(tunney$tr_gene) )
nottraining = which( !as.logical(tunney$tr_gene) )

liu.subspace = read.delim( liu_s_fname, header=F)
liu.names = read.delim( liu_gene_fname, header=F)
row.names(liu.subspace) = liu.names$V1

#### combine them all -- genes that are found in all
in.all = intersect( tunney.names[nottraining], row.names(liu.subspace) )

ord = order( tunney[in.all,"fp_density"], decreasing=T )
in.all = in.all[ord] # sort gene names by density

num.all = length(in.all)

tunney.x = match( in.all, row.names(tunney) )
liu.sub.x = match( in.all, row.names(liu.subspace) )

# tunney.loess = predict( loess( tunney$pearson_r[tunney.x] ~ c(1:num.all) ),
#                         1:num.all, se=T )
# liu.sub.loess = lapply( 1:8, function(i){predict( loess( liu.subspace[liu.sub.x,i] ~ c(1:num.all) ),
#                          1:num.all, se=T )} ) 


cnames = c(paste0( "Liu V", 0:7), "Tunney")

all.data = as.data.frame( cbind(liu.subspace[liu.sub.x,],
                                tunney$pearson_r[tunney.x] )
)
names( all.data ) = cnames

# loess.fits = as.data.frame( cbind( sapply( liu.sub.loess, function(i) {i$fit}),
#                                    tunney.loess$fit ),
#                             col.names = cnames
# )

cols = c( rev(gray.colors(8)), "red")

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

pdf( out_fname, width=3, height=3, pointsize=7, useDingbats=F )
par( mex = 0.65 ) # sets margin stuff
par( mar = c(6,5.5,5,3) )
par( oma = c(0,1.5,1,0) )
par( lwd = 0.75 )
par( xpd = NA )

plot( NA,
      xlim = c(1, num.all),
      ylim = c(0,1),
      bty = "n",
#      axes = F,
      xlab = "genes sorted by footprint density",
      ylab = "Pearson correlation per gene"
)
# plot all loess lines
Map( function(x,y){ lines( 1:num.all, x, col=y) }, loess.fits, cols )

#axis( 1 )
#axis( 2, at = c(0,0.5,1), labels = c(0, NA, 1), las = 1)

legend( "topright", inset=c( 0, 0 ), 
        legend = cnames[5:9],
        fill = cols[5:9],
        border = NA,
        bty="n")

legend( "topright", inset=c( 0.3, 0), 
        legend = cnames[1:4],
        fill = cols[1:4],
        border = NA,
        bty="n")

dev.off()
