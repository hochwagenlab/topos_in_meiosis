# Topo figure S2
####################################################################################
####################################################################################
# Figure S2a

Top1_myc = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/AH9847Myc-3h-735-841-Reps-SK1Yue-B3W4-MACS2/AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top1_mycd = gendiv(Top1_myc)

gff <- hwglabr2::get_gff('SK1Yue')
transcription <- read.csv('/Volumes/LabShare/HTGenomics/HiSeqOutputs/RNA-seq/2016.03.16-2h+3h/2017.06.16_SK1Yue_EdgeR_tpm.csv')
gff <- gff[which(gff$type=='gene')]
colnames(transcription)[1] <- "ID"
gff <- data.frame(gff)
gff_txn <- merge(x=gff,y=transcription[,c(1,3)],by='ID', all.x = TRUE)
gff_txn <- gff_txn[which(gff_txn$seqnames!='chrMT'),]
gff_txn <- gff_txn[which(gff_txn$seqnames!='scplasm1'),]
gffgr <- GRanges(seqnames = gff_txn$seqnames,IRanges(start = gff_txn$start,end=gff_txn$end),strand = gff_txn$strand,transcription=gff_txn$AH119_3h,gene=gff_txn$Name)
negstrand=gffgr[strand(gffgr)=='-']
posstrand=gffgr[strand(gffgr)=='+']
negstrandswitch <- GRanges(seqnames = seqnames(negstrand),ranges = IRanges(start=end(negstrand),end=end(negstrand)),strand=strand(negstrand),transcription=negstrand$transcription,gene=negstrand$gene)
end(posstrand) <- start(posstrand)
gffall <- c(posstrand,negstrandswitch)
gffall_sort <- gffall[order(gffall$transcription,decreasing = T)]
gffall_sort <- gffall_sort[which(gffall_sort$transcription!='<NA>')]

Top1_atg <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, gffall_sort,
                                               extend=c(500,0), w=1,empty_value=NA,
                                               mean_mode="weighted",
                                               value_column='score')
par(las=1)
sigrna <- data.frame(cbind(gffall_sort$transcription,rowSums(data.frame(Top1_atg))))
colnames(sigrna) <- c('txn','top1')
plot(log(sigrna$txn),log(sigrna$top1))
abline(lm(data=sigrna,log(top1)~log(txn)))
summary(lm(data=sigrna,log(top1)~log(txn)))
cor(x=sigrna$txn,y=sigrna$top1, use="complete.obs", method="pearson")
#[1] 0.2055718
####################################################################################
# Figure S2b

Top2_wt = hwglabr2::import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-413-504-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top2_wtd = gendiv(Top2_wt)

gff <- hwglabr2::get_gff('SK1Yue')
transcription <- read.csv('/Volumes/LabShare/HTGenomics/HiSeqOutputs/RNA-seq/2016.03.16-2h+3h/2017.06.16_SK1Yue_EdgeR_tpm.csv')
gff <- gff[which(gff$type=='gene')]
colnames(transcription)[1] <- "ID"
gff <- data.frame(gff)
gff_txn <- merge(x=gff,y=transcription[,c(1,3)],by='ID', all.x = TRUE)
gff_txn <- gff_txn[which(gff_txn$seqnames!='chrMT'),]
gff_txn <- gff_txn[which(gff_txn$seqnames!='scplasm1'),]
gffgr <- GRanges(seqnames = gff_txn$seqnames,IRanges(start = gff_txn$start,end=gff_txn$end),strand = gff_txn$strand,transcription=gff_txn$AH119_3h,gene=gff_txn$Name)
negstrand=gffgr[strand(gffgr)=='-']
posstrand=gffgr[strand(gffgr)=='+']

negstrandswitch <- GRanges(seqnames = seqnames(negstrand),ranges = IRanges(start=end(negstrand),end=end(negstrand)),strand=strand(negstrand),transcription=negstrand$transcription,gene=negstrand$gene)
end(posstrand) <- start(posstrand)
gffall <- c(posstrand,negstrandswitch)
gffall_sort <- gffall[order(gffall$transcription,decreasing = T)]
gffall_sort <- gffall_sort[which(gffall_sort$transcription!='<NA>')]
Top2_atg <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, gffall_sort,
                                               extend=c(500,0), w=10,empty_value=NA,
                                               mean_mode="weighted",
                                               value_column='score')

par(las=1)
sigrna2 <- data.frame(cbind(gffall_sort$transcription,rowSums(data.frame(Top2_atg))))
colnames(sigrna2) <- c('txn','top2')
plot(log(sigrna2$txn),log(sigrna2$top2))
abline(lm(data=sigrna2,log(top2)~log(txn)))
summary(lm(data=sigrna2,log(top2)~log(txn)))
cor(x=sigrna2$txn,y=sigrna2$top2, use="complete.obs", method="pearson")
#[1] 0.1829016

####################################################################################
# Figure S2c
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
mean(lowlow$width) #936.1322
mean(highhigh$width) #929.5546
median(lowlow$width) #434
median(highhigh$width) #733

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

mean(lowlow$width) #795.7273
mean(highhigh$width) #1074.407
median(lowlow$width) #407
median(highhigh$width) #772.5

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

mean(lowlow$width) #1043.414
mean(highhigh$width) #845.4915
median(lowlow$width) #478
median(highhigh$width) #720
