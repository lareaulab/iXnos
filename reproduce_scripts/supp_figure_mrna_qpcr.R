args <- commandArgs(trailingOnly = TRUE)
citrine_construct_scores_fname = args[1]
mrna_data_fname1 = args[2]
mrna_data_fname2 = args[3]
mrna_data_fname3 = args[4]
out_file = args[5]

cit = read.delim(citrine_construct_scores_fname, header=T, row.names = 1, comment.char="#")

qpcr.data = read.delim(mrna_data_fname1, sep=",", header=T )
qpcr.data = qpcr.data[which(!(qpcr.data$Strain == "MAX" & qpcr.data$Isolate == 3)),]
## this sample's pellet was aspirated - removing from all analysis

qpcr.2.data = read.delim(mrna_data_fname2, sep=",", header=T )
qpcr.2.data = subset( qpcr.2.data, strain != "Y66")
## 666 didn't grow overnight in this experiment - drop that strain

qpcr.3.data = read.delim(mrna_data_fname3, sep=",", header=T )

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

## paired mcherry and citrine amplifications. roughly, ratio of starting citrine:mcherry amounts
mrna.ratio = data.frame(mch.amp / cit.amp)
## ADD IN A NORMALIZATION TO THE MIN MCHERRY:CITRINE RATIO
mrna.ratio = mrna.ratio / median( mrna.ratio[,"MIN"])

## second experiment (min, 0, 6, 9)
cit.2.cy0 = t( tapply(qpcr.2.data$Cy0, list(qpcr.2.data$strain, qpcr.2.data$isolate, qpcr.2.data$amplicon2), median)[,,"eCitrine"] )
mch.2.cy0 = t( tapply(qpcr.2.data$Cy0, list(qpcr.2.data$strain, qpcr.2.data$isolate, qpcr.2.data$amplicon2), median)[,,"mCherry"] )
cit.2.effs = tapply( qpcr.2.data$eff, list(qpcr.2.data$strain,  qpcr.2.data$amplicon2), median)[,"eCitrine"]
mch.2.eff = tapply( qpcr.2.data$eff, qpcr.2.data$amplicon2, median)["mCherry"]
mch.2.amp = mch.2.eff ^ mch.2.cy0
cit.2.amp = sapply(colnames(cit.2.cy0), function(x){ cit.2.effs[x] ^ cit.2.cy0[,x]})
mrna.2.ratio = data.frame(mch.2.amp / cit.2.amp)
mrna.2.ratio = mrna.2.ratio / median( mrna.2.ratio[,"MIN"])
mrna.2.ratio = subset( mrna.2.ratio, select = -Y66)
colnames(mrna.2.ratio) = c("MIN","Y000","Y999")
## remove Y666 because it didn't grow in this experiment

## third experiment (min, 6)
cit.3.cy0 = t( tapply(qpcr.3.data$Cy0, list(qpcr.3.data$strain, qpcr.3.data$isolate, qpcr.3.data$amplicon2), median)[,,"eCitrine"] )
mch.3.cy0 = t( tapply(qpcr.3.data$Cy0, list(qpcr.3.data$strain, qpcr.3.data$isolate, qpcr.3.data$amplicon2), median)[,,"mCherry"] )
cit.3.effs = tapply( qpcr.3.data$eff, list(qpcr.3.data$strain,  qpcr.3.data$amplicon2), median)[,"eCitrine"]
mch.3.eff = tapply( qpcr.3.data$eff, qpcr.3.data$amplicon2, median)["mCherry"]
mch.3.amp = mch.3.eff ^ mch.3.cy0
cit.3.amp = sapply(colnames(cit.3.cy0), function(x){ cit.3.effs[x] ^ cit.3.cy0[,x]})
mrna.3.ratio = data.frame(mch.3.amp / cit.3.amp)
mrna.3.ratio = mrna.3.ratio / median( mrna.3.ratio[,"MIN"])
colnames(mrna.3.ratio) = c("MIN","Y666")

## combine experiments (normalized to MIN in each experiment)
## only keep the MIN measurements from the first experiment
mrna.ratios = cbind( mrna.ratio, subset(mrna.2.ratio, select=-MIN), subset(mrna.3.ratio, select=-MIN))

nn.scores =  sapply( colnames(mrna.ratios), function(x) { rep( cit[x,"nn.score"], 3 )})

collist = list( MIN="magenta3", CHA2="purple2", Y000="royalblue2", Y333="green3", Y666="gold1", Y999="darkorange2", MAX="red2")
cols = rep( unlist(collist[colnames(mrna.ratios)]), each=3 )

pdf(out_file, width=2, height=1.67, pointsize=7, useDingbats=F, bg="white" )
#cairo_pdf(out_file, width=2, height=1.67, pointsize=7 )
par( mex = 0.65 ) # sets margin stuff
par( mar =c(7,6.5,4,3) )
par( oma = c(0,0.5,1,0) )
plot(as.numeric(nn.scores), unlist(mrna.ratios),
     col = cols,
     cex = 0.4,
     axes = F,
     xlim = c(150,400),
     ylim = c(0, max(mrna.ratio, na.rm=T)),
#     pch = 3,
     pch = 20,
     xlab = "",
     ylab = "eCitrine / mCherry\nrelative mRNA ratio"
)
axis( 1 )
axis( 2 )
title( xlab = "predicted elongation time\n(arbitrary units)", line = 4.5 )
mtext( "B", font = 2, line = -3, side = 3, outer = T, adj = 0 ) 
dev.off()
