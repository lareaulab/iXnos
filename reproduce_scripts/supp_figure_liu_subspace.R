# analysis of Liu subspace performance
# Liu data from running riboshape on the Weinberg dataset locally

args <- commandArgs(trailingOnly = TRUE)
tunney_fname = args[1] # corrs_by_gene_density.txt
liu_s_fname = args[2] # ASK: subspace_corrs.txt
liu_gene_fname = args[3] # ASK: gene_names.txt
out_fname = args[4] 

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

all.data = as.data.frame( cbind(liu.subspace[liu.sub.x,],
                                tunney$pearson_r[tunney.x] )
)
names( all.data ) = c(paste0( "LiuV", 0:7), "Tunney")

# pdf( out_fname, width=4, height=3, pointsize=7, useDingbats=F )
# par( mex = 0.65 ) # sets margin stuff
# par( mar = c(6,5.5,5,3) )
# par( oma = c(0,1.5,1,0) )
# par( lwd = 0.75 )
# par( xpd = F )
# boxplot( all.data,
#          notch = T,
#          bty = "n",
#          ylim = c(-1,1),
#          frame = F,
#          col = "grey",
#          boxwex = 0.6,
#          whisklty = 1,
#          axes = F
# )
# axis(2)
# axis(1, at = 1:9, 
#      lwd = 0,
#      labels = c( expression("Liu V"[0]), 
#                     expression("V"[1]),
#                     expression("V"[2]),
#                     expression("V"[3]),
#                     expression("V"[4]),
#                     expression("V"[5]),
#                     expression("V"[6]),
#                     expression("V"[7]),
#                     "Tunney"))
# abline( h = median( all.data$Tunney), lty = 2, col="red")
# 
# dev.off()

round( sapply(all.data, mean), 2 )
