args <- commandArgs(trailingOnly = TRUE)
citrine_construct_scores_fname = args[1]
mrna_data_fname = args[2]
out_fname = args[3]

## mRNA
cit = read.delim(citrine_construct_scores_fname, header=T, row.names = 1, comment.char="#")

qpcr.data = read.delim(mrna_data_fname, sep=",", header=T )
qpcr.data = qpcr.data[which(!(qpcr.data$Strain == "MAX" & qpcr.data$Isolate == 3)),]
# this sample's pellet was aspirated - removing from all analysis

# median cy0 for each group of three technical replicates (so each strain+isolate gets one value for cit and one for mch)
cit.cy0 = t( tapply(qpcr.data$Cy0, list(qpcr.data$Strain, qpcr.data$Isolate, qpcr.data$amplicon2), median)[,,"eCitrine"] )
mch.cy0 = t( tapply(qpcr.data$Cy0, list(qpcr.data$Strain, qpcr.data$Isolate, qpcr.data$amplicon2), median)[,,"mCherry"] )

# median efficiency across all times we used that citrine primer pair (3 isolates x 3 tech reps)
cit.effs = tapply( qpcr.data$eff, list(qpcr.data$Strain,  qpcr.data$amplicon2), median)[,"eCitrine"]
# median efficiency across all mCherry qPCRs (4 strains x 3 isolates x 3 tech reps)
mch.eff = tapply( qpcr.data$eff, qpcr.data$amplicon2, median)["mCherry"]

# mcherry amplification for each strain+isolate: median efficiency for all mch experiments, median mcherry cy0 for each strain+isolate
mch.amp = mch.eff ^ mch.cy0
# citrine amplification for each strain+isolate: median efficiency for each cit primer set, median citrine cy0 for each strain+isolate
cit.amp = sapply(colnames(cit.cy0), function(x){ cit.effs[x] ^ cit.cy0[,x]})

# paired mcherry and citrine amplifications. roughly, ratio of starting citrine:mcherry amounts
mrna.ratio = mch.amp / cit.amp
# ADD IN A NORMALIZATION TO THE CHA MCHERRY:CITRINE RATIO
mrna.ratio = mrna.ratio / median( mrna.ratio[,"CHA2"])

nn.scores =  sapply( colnames(mrna.ratio), function(x) { rep( cit$nn.score[x], 3 )}) 

cols = rep( c("magenta","red","purple","cyan"), each = 3)

# protein
facs.data = read.delim("~/Berkeley/RegressionPaperLL/CodOpt/remeasuring_diploids_20170829/gated-facs-data.csv", sep=",", header=T )
# median green/red ratio of all cells included in analysis 
facs.medians = data.frame(tapply( facs.data$green / facs.data$red, list(facs.data$Isolate, facs.data$Strain), median ))
facs.medians = subset(facs.medians, select = c(CHA2, MAX, MIN, Y333))
facs.medians = facs.medians[c(3,5,6),]

# one per strain+isolate
te = facs.medians/mrna.ratio

cairo_pdf("figure_3D_te.pdf", width=2, height=1.67, pointsize=7 )
par( mex = 0.65 ) # sets margin stuff
par( mar =c(6,6.5,5,3) )
par( oma = c(0,0.5,1,0) )
par( xpd = NA )
plot(as.numeric(nn.scores), unlist(te),
     col = cols,
     cex = 0.6,
     axes = F,
     xlim = c(150,400),
     ylim = c(0, max(te,na.rm=T)),
     pch = 3,
#     pch = 20,
     xlab = "NN score",
     ylab = "translation efficiency"
)
axis( 1 )
axis( 2 )
#axis( 2, at = seq( 0, 0.8, by = 0.2 ))
mtext( "D", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 
dev.off()
