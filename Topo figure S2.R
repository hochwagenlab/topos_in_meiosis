# Topo figure S2
####################################################################################
####################################################################################
# Figure S2
library(hwglabr2)
library(GenomicRanges)

brartotmrna <- read.csv('/Volumes/LabShare/Jonna/papers/Topo/figures/RNAseq/GSE108778_timecourse_replicate_2_totRNA.txt.gz',sep='\t',header=T)
mrna3h <- brartotmrna[,c('gene','X3hr.totRNA.rpkm')]
mrna10h <- brartotmrna[,c('gene','X10hr.totRNA.rpkm')]
mrnaAA <- brartotmrna[,c('gene','MATa.a.totRNA.rpkm')]

gff <- hwglabr2::get_gff('SK1Yue')
intergen <- hwglabr2::get_intergenic_regions('SK1Yue',as_gr = F)
divergent <- intergen[which(intergen$type=='divergent'),]

# Meiotic 3h

txn_data <- 'X3hr.totRNA.rpkm'
transcription <- mrna3h
gff <- gff[which(gff$type=='gene')]
colnames(transcription)[1] <- "Name"
gff <- data.frame(gff)
gff_txn <- merge(x=gff,y=transcription,by='Name', all.x = TRUE)
gff_txn <- gff_txn[which(gff_txn$seqnames!='chrMT'),]
gff_txn <- gff_txn[which(gff_txn$seqnames!='scplasm1'),]

colnames(gff_txn)[1] <- 'left_gene'
divergent_txn <- merge(x=divergent,y=gff_txn[,c('left_gene',txn_data)],by='left_gene', all.x = TRUE)
colnames(divergent_txn)[10] <- 'left_gene_txn'
colnames(gff_txn)[1] <- 'right_gene'
divergent_txn <- merge(x=divergent_txn,y=gff_txn[,c('right_gene',txn_data)],by='right_gene', all.x = TRUE)
colnames(divergent_txn)[11] <- 'right_gene_txn'

divergent_txn_sort <- divergent_txn[order(divergent_txn$left_gene_txn),]
highhigh <- divergent_txn_sort[which(divergent_txn_sort$left_gene_txn>=70&divergent_txn_sort$right_gene_txn>=70),]
lowlow <- divergent_txn_sort[which(divergent_txn_sort$left_gene_txn<10&divergent_txn_sort$right_gene_txn<10),]
nrow(highhigh);nrow(lowlow)

highhigh[,'width'] <- highhigh$right_coordinate-highhigh$left_coordinate+1
lowlow[,'width'] <- lowlow$right_coordinate-lowlow$left_coordinate+1
alldata <- list(lowlow$width,highhigh$width)
par(las=1)
boxplot(alldata,ylab="Intergenic size (bp)",names=c('low txn','high txn'),frame.plot=F,cex.lab=1.5,cex.axis=1.25,outline=F)
wilcox.test(lowlow$width,highhigh$width,paired=F) # p-value = 1.131e-05

# Meiotic 10h

txn_data <- 'X10hr.totRNA.rpkm'
transcription <- mrna10h
gff <- gff[which(gff$type=='gene')]
colnames(transcription)[1] <- "Name"
gff <- data.frame(gff)
gff_txn <- merge(x=gff,y=transcription,by='Name', all.x = TRUE)
gff_txn <- gff_txn[which(gff_txn$seqnames!='chrMT'),]
gff_txn <- gff_txn[which(gff_txn$seqnames!='scplasm1'),]

colnames(gff_txn)[1] <- 'left_gene'
divergent_txn <- merge(x=divergent,y=gff_txn[,c('left_gene',txn_data)],by='left_gene', all.x = TRUE)
colnames(divergent_txn)[10] <- 'left_gene_txn'
colnames(gff_txn)[1] <- 'right_gene'
divergent_txn <- merge(x=divergent_txn,y=gff_txn[,c('right_gene',txn_data)],by='right_gene', all.x = TRUE)
colnames(divergent_txn)[11] <- 'right_gene_txn'

divergent_txn_sort <- divergent_txn[order(divergent_txn$left_gene_txn),]
highhigh <- divergent_txn_sort[which(divergent_txn_sort$left_gene_txn>=60&divergent_txn_sort$right_gene_txn>=60),]
lowlow <- divergent_txn_sort[which(divergent_txn_sort$left_gene_txn<10&divergent_txn_sort$right_gene_txn<10),]
nrow(highhigh);nrow(lowlow)

highhigh[,'width'] <- highhigh$right_coordinate-highhigh$left_coordinate+1
lowlow[,'width'] <- lowlow$right_coordinate-lowlow$left_coordinate+1
alldata <- list(lowlow$width,highhigh$width)
par(las=1)
boxplot(alldata,ylab="Intergenic size (bp)",names=c('low txn','high txn'),frame.plot=F,cex.lab=1.5,cex.axis=1.25,outline=F)
wilcox.test(lowlow$width,highhigh$width,paired=F) # p-value = 3.933e-11

# Non-meiotic

txn_data <- 'MATa.a.totRNA.rpkm'
transcription <- mrnaAA
gff <- gff[which(gff$type=='gene')]
colnames(transcription)[1] <- "Name"
gff <- data.frame(gff)
gff_txn <- merge(x=gff,y=transcription,by='Name', all.x = TRUE)
gff_txn <- gff_txn[which(gff_txn$seqnames!='chrMT'),]
gff_txn <- gff_txn[which(gff_txn$seqnames!='scplasm1'),]

colnames(gff_txn)[1] <- 'left_gene'
divergent_txn <- merge(x=divergent,y=gff_txn[,c('left_gene',txn_data)],by='left_gene', all.x = TRUE)
colnames(divergent_txn)[10] <- 'left_gene_txn'
colnames(gff_txn)[1] <- 'right_gene'
divergent_txn <- merge(x=divergent_txn,y=gff_txn[,c('right_gene',txn_data)],by='right_gene', all.x = TRUE)
colnames(divergent_txn)[11] <- 'right_gene_txn'

divergent_txn_sort <- divergent_txn[order(divergent_txn$left_gene_txn),]
highhigh <- divergent_txn_sort[which(divergent_txn_sort$left_gene_txn>=60&divergent_txn_sort$right_gene_txn>=60),]
lowlow <- divergent_txn_sort[which(divergent_txn_sort$left_gene_txn<9&divergent_txn_sort$right_gene_txn<9),]
nrow(highhigh);nrow(lowlow)

highhigh[,'width'] <- highhigh$right_coordinate-highhigh$left_coordinate+1
lowlow[,'width'] <- lowlow$right_coordinate-lowlow$left_coordinate+1
alldata <- list(lowlow$width,highhigh$width)
par(las=1)
boxplot(alldata,ylab="Intergenic size (bp)",names=c('low txn','high txn'),frame.plot=F,cex.lab=1.5,cex.axis=1.25,outline=F)
wilcox.test(lowlow$width,highhigh$width,paired=F) # p-value = 0.004679