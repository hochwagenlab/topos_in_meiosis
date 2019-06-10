# Topo figure 2
####################################################################################
####################################################################################
library(EnrichedHeatmap)
library(circlize)
library(GenomicRanges)
library(hwglabr2)
library(tidyverse)

# Figure 2a

Top1_myc = import_bedGraph("AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top1_mycd = gendiv(Top1_myc)

gff <- hwglabr2::get_gff('SK1Yue')
transcription <- read.csv('2017.06.16_SK1Yue_EdgeR_tpm.csv')
gff <- gff[which(gff$type=='gene')]
colnames(transcription)[1] <- "ID"
gff <- data.frame(gff)
gff_txn <- merge(x=gff,y=transcription[,c(1,3)],by='ID', all.x = TRUE)
gff_txn <- gff_txn[which(gff_txn$seqnames!='chrMT'),]
gff_txn <- gff_txn[which(gff_txn$seqnames!='scplasm1'),]

# ATG heatmap for Top1

gffgr <- GRanges(seqnames = gff_txn$seqnames,IRanges(start = gff_txn$start,end=gff_txn$end),strand = gff_txn$strand,transcription=gff_txn$AH119_3h,gene=gff_txn$Name)
negstrand=gffgr[strand(gffgr)=='-']
posstrand=gffgr[strand(gffgr)=='+']
negstrandswitch <- GRanges(seqnames = seqnames(negstrand),ranges = IRanges(start=end(negstrand),end=end(negstrand)),strand=strand(negstrand),transcription=negstrand$transcription,gene=negstrand$gene)
end(posstrand) <- start(posstrand)
gffall <- c(posstrand,negstrandswitch)
gffall_sort <- gffall[order(gffall$transcription,decreasing = T)]
gffall_sort <- gffall_sort[which(gffall_sort$transcription!='<NA>')]
mcols(gffall_sort) <- DataFrame(class=c(rep(1:4, each=length(gffall_sort)/4),4,4))

Top1_atg <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, gffall_sort,
                                               extend=c(500,500), w=10,empty_value=NA,
                                               mean_mode="weighted",
                                               value_column='score')
col_fun <- colorRamp2(quantile(Top1_atg, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
partition <- gffall_sort$class
EnrichedHeatmap(Top1_atg, col = col_fun, name = "Top1", row_title_rot = 0,
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:4),
                                                                         show_error = T)),
                top_annotation_height = unit(2, "cm"),
                row_order = 1:length(gffall_sort),
                split=gffall_sort$class,
                axis_name = c("-500", "START","500"))+
  Heatmap(partition, col = structure(1:4, names = as.character(1:4)), name = "",row_order = 1:length(gffall_sort),
          show_row_names = FALSE, width = unit(5, "mm"))

## STOP heatmap for Top2

gffgr <- GRanges(seqnames = gff_txn$seqnames,IRanges(start = gff_txn$start,end=gff_txn$end),strand = gff_txn$strand,transcription=gff_txn$AH119_3h,gene=gff_txn$Name)
negstrand=gffgr[strand(gffgr)=='-']
posstrand=gffgr[strand(gffgr)=='+']
negstrandswitch <- GRanges(seqnames = seqnames(negstrand),ranges = IRanges(start=start(negstrand),end=start(negstrand)),strand=strand(negstrand),transcription=negstrand$transcription,gene=negstrand$gene)
start(posstrand) <- end(posstrand)
gffall <- c(posstrand,negstrandswitch)
gffall_sort <- gffall[order(gffall$transcription,decreasing = T)]
gffall_sort <- gffall_sort[which(gffall_sort$transcription!='<NA>')]
mcols(gffall_sort) <- DataFrame(class=c(rep(1:4, each=length(gffall_sort)/4),4,4))

Top1_stop <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, gffall_sort,
                                                extend=c(500,500), w=10,empty_value=NA,
                                                mean_mode="weighted",
                                                value_column='score')
col_fun <- colorRamp2(quantile(Top1_atg, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
partition <- gffall_sort$class
EnrichedHeatmap(Top1_stop, col = col_fun, name = "Top1", row_title_rot = 0,
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:4),
                                                                         show_error = T)),
                top_annotation_height = unit(2, "cm"),
                row_order = 1:length(gffall_sort),
                split=gffall_sort$class,
                axis_name = c("-500", "STOP","500"))+
  Heatmap(partition, col = structure(1:4, names = as.character(1:4)), name = "",row_order = 1:length(gffall_sort),
          show_row_names = FALSE, width = unit(5, "mm"))
####################################################################################
####################################################################################
# Figure 2b

Top2_wt = hwglabr2::import_bedGraph("Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top2_wtd = gendiv(Top2_wt)

gff <- hwglabr2::get_gff('SK1Yue')
transcription <- read.csv('2017.06.16_SK1Yue_EdgeR_tpm.csv')
gff <- gff[which(gff$type=='gene')]
colnames(transcription)[1] <- "ID"
gff <- data.frame(gff)
gff_txn <- merge(x=gff,y=transcription[,c(1,3)],by='ID', all.x = TRUE)
gff_txn <- gff_txn[which(gff_txn$seqnames!='chrMT'),]
gff_txn <- gff_txn[which(gff_txn$seqnames!='scplasm1'),]

# ATG heatmap for Top2
gffgr <- GRanges(seqnames = gff_txn$seqnames,IRanges(start = gff_txn$start,end=gff_txn$end),strand = gff_txn$strand,transcription=gff_txn$AH119_3h,gene=gff_txn$Name)
negstrand=gffgr[strand(gffgr)=='-']
posstrand=gffgr[strand(gffgr)=='+']

negstrandswitch <- GRanges(seqnames = seqnames(negstrand),ranges = IRanges(start=end(negstrand),end=end(negstrand)),strand=strand(negstrand),transcription=negstrand$transcription,gene=negstrand$gene)
end(posstrand) <- start(posstrand)
gffall <- c(posstrand,negstrandswitch)
gffall_sort <- gffall[order(gffall$transcription,decreasing = T)]
gffall_sort <- gffall_sort[which(gffall_sort$transcription!='<NA>')]
mcols(gffall_sort) <- DataFrame(class=c(rep(1:4, each=length(gffall_sort)/4),4,4))

Top2_atg <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, gffall_sort,
                                               extend=c(500,500), w=10,empty_value=NA,
                                               mean_mode="weighted",
                                               value_column='score')
col_fun <- colorRamp2(quantile(Top2_atg, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
partition <- gffall_sort$class
EnrichedHeatmap(Top2_atg, col = col_fun, name = "Top2", row_title_rot = 0,
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:4),
                                                                         show_error = T)),
                top_annotation_height = unit(2, "cm"),
                row_order = 1:length(gffall_sort),
                split=gffall_sort$class,
                axis_name = c("-500", "START","500"))+
  Heatmap(partition, col = structure(1:4, names = as.character(1:4)), name = "",row_order = 1:length(gffall_sort),
          show_row_names = FALSE, width = unit(5, "mm"))

## STOP heatmap for Top2

gffgr <- GRanges(seqnames = gff_txn$seqnames,IRanges(start = gff_txn$start,end=gff_txn$end),strand = gff_txn$strand,transcription=gff_txn$AH119_3h,gene=gff_txn$Name)
negstrand=gffgr[strand(gffgr)=='-']
posstrand=gffgr[strand(gffgr)=='+']
negstrandswitch <- GRanges(seqnames = seqnames(negstrand),ranges = IRanges(start=start(negstrand),end=start(negstrand)),strand=strand(negstrand),transcription=negstrand$transcription,gene=negstrand$gene)
start(posstrand) <- end(posstrand)
gffall <- c(posstrand,negstrandswitch)
gffall_sort <- gffall[order(gffall$transcription,decreasing = T)]
gffall_sort <- gffall_sort[which(gffall_sort$transcription!='<NA>')]
mcols(gffall_sort) <- DataFrame(class=c(rep(1:4, each=length(gffall_sort)/4),4,4))

Top2_stop <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, gffall_sort,
                                                extend=c(500,500), w=10,empty_value=NA,
                                                mean_mode="weighted",
                                                value_column='score')
col_fun <- colorRamp2(quantile(Top2_atg, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
partition <- gffall_sort$class
EnrichedHeatmap(Top2_stop, col = col_fun, name = "Top2", row_title_rot = 0,
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:4),
                                                                         show_error = T)),
                top_annotation_height = unit(2, "cm"),
                row_order = 1:length(gffall_sort),
                split=gffall_sort$class,
                axis_name = c("-500", "STOP","500"))+
  Heatmap(partition, col = structure(1:4, names = as.character(1:4)), name = "",row_order = 1:length(gffall_sort),
          show_row_names = FALSE, width = unit(5, "mm"))

####################################################################################
####################################################################################
# Figure 2c

# Divergent IGRs
gff <- hwglabr2::get_gff('SK1Yue')
transcription <- read.csv('2017.06.16_SK1Yue_EdgeR_tpm.csv')
intergen <- hwglabr2::get_intergenic_regions('SK1Yue',as_gr = T)
divergent <- intergen[which(intergen$type=='divergent')]
mcols(divergent)['widths'] <- width(divergent)
divergent <- data.frame(divergent)
gff <- gff[which(gff$type=='gene')]
colnames(transcription)[1] <- "ID"
gff <- data.frame(gff)
gff_txn <- merge(x=gff,y=transcription[,c(1,3)],by='ID', all.x = TRUE)
gff_txn <- gff_txn[which(gff_txn$seqnames!='chrMT'),]
gff_txn <- gff_txn[which(gff_txn$seqnames!='scplasm1'),]

colnames(gff_txn)[11] <- 'left_gene'
divergent_txn <- merge(x=divergent,y=gff_txn[,c('left_gene','AH119_3h')],by='left_gene', all.x = TRUE)
colnames(divergent_txn)[12] <- 'left_gene_txn'
colnames(gff_txn)[11] <- 'right_gene'
divergent_txn <- merge(x=divergent_txn,y=gff_txn[,c('right_gene','AH119_3h')],by='right_gene', all.x = TRUE)
colnames(divergent_txn)[13] <- 'right_gene_txn'
divergent_txn_uniq <- divergent_txn[!duplicated(divergent_txn$left_gene),]
divergent_txn_uniq <- divergent_txn_uniq[!duplicated(divergent_txn_uniq$right_gene),]
divergent_txn_sort <- divergent_txn_uniq[order(divergent_txn_uniq$left_gene_txn),]
highhigh <- divergent_txn_sort[which(divergent_txn_sort$left_gene_txn>=125&divergent_txn_sort$right_gene_txn>=125),]
lowlow <- divergent_txn_sort[which(divergent_txn_sort$left_gene_txn<50&divergent_txn_sort$right_gene_txn<50),]
nrow(highhigh);nrow(lowlow)

allwidths <- list(lowlow$widths,highhigh$widths)
boxplot(allwidths,names=c('low txn','high txn'),ylab="Gene width (bp)",frame.plot=F,cex.lab=1.5,cex.axis=1.25,outline=F)
wilcox.test(highhigh$widths,lowlow$widths,paired=F) # p-value = 2.736e-08***

# Tandem IGRs

gff <- hwglabr2::get_gff('SK1Yue')
transcription <- read.csv('2017.06.16_SK1Yue_EdgeR_tpm.csv')
intergen <- hwglabr2::get_intergenic_regions('SK1Yue',as_gr = T)
tandem <- intergen[which(intergen$type=='tandem')]
mcols(tandem)['widths'] <- width(tandem)
tandem <- data.frame(tandem)

gff <- gff[which(gff$type=='gene')]
colnames(transcription)[1] <- "ID"
gff <- data.frame(gff)
gff_txn <- merge(x=gff,y=transcription[,c(1,3)],by='ID', all.x = TRUE)
gff_txn <- gff_txn[which(gff_txn$seqnames!='chrMT'),]
gff_txn <- gff_txn[which(gff_txn$seqnames!='scplasm1'),]

colnames(gff_txn)[11] <- 'left_gene'
tandem_txn <- merge(x=tandem,y=gff_txn[,c('left_gene','AH119_3h')],by='left_gene', all.x = TRUE)
colnames(tandem_txn)[12] <- 'left_gene_txn'
colnames(gff_txn)[11] <- 'right_gene'
tandem_txn <- merge(x=tandem_txn,y=gff_txn[,c('right_gene','AH119_3h')],by='right_gene', all.x = TRUE)
colnames(tandem_txn)[13] <- 'right_gene_txn'
tandem_txn <- tandem_txn[which(tandem_txn$left_gene_txn!='NA'),]
tandem_txn <- tandem_txn[which(tandem_txn$right_gene_txn!='NA'),]
tandem_txn_plus <- tandem_txn[which(tandem_txn$left_gene_strand=='+'),]
colnames(tandem_txn_plus) <- c('right_gene','left_gene','chr','end','start','width','strand','type','usstrand','dsstrand','widths','left_gene_txn','right_gene_txn')
tandem_txn_minus <- tandem_txn[which(tandem_txn$left_gene_strand=='-'),]
colnames(tandem_txn_minus) <- c('right_gene','left_gene','chr','end','start','width','strand','type','usstrand','dsstrand','widths','right_gene_txn','left_gene_txn')
tandem_txn_sign <- rbind(tandem_txn_plus,tandem_txn_minus)

tandem_txn_uniq <- tandem_txn_sign[!duplicated(tandem_txn_sign$left_gene),]
tandem_txn_uniq <- tandem_txn_uniq[!duplicated(tandem_txn_uniq$right_gene),]
highhigh <- tandem_txn_sign[which(tandem_txn_sign$left_gene_txn>=150&tandem_txn_sign$right_gene_txn>=150),] #saved images: 125 &50
lowlow <- tandem_txn_sign[which(tandem_txn_sign$left_gene_txn<40&tandem_txn_sign$right_gene_txn<40),]
nrow(highhigh);nrow(lowlow)

allwidths <- list(lowlow$widths,highhigh$widths)
boxplot(allwidths,names=c('low txn','high txn'),ylab="Gene width (bp)",frame.plot=F,cex.lab=1.5,cex.axis=1.25,outline=F)
wilcox.test(highhigh$widths,lowlow$widths,paired=F) # p-value = 0.7127

####################################################################################
####################################################################################
# Figure 2g
brartotmrna <- read.csv('GSE108778_timecourse_replicate_2_totRNA.txt.gz',sep='\t',header=T)
mrna0h <- brartotmrna[,c('gene','X0hr.totRNA.rpkm')]
mrnaveg <- brartotmrna[,c('gene','vexp.repl1.totrna.rpkm')]

gff <- hwglabr2::get_gff('SK1Yue')
gff <- gff[which(gff$type=='gene')]
intergen <- hwglabr2::get_intergenic_regions('SK1Yue',as_gr = F)
divergent <- intergen[which(intergen$type=='divergent'),]

# Vegetative

txn_data <- 'vexp.repl1.totrna.rpkm'
transcription <- mrnaveg
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
highhigh <- divergent_txn_sort[which(divergent_txn_sort$left_gene_txn>=230&divergent_txn_sort$right_gene_txn>=230),]
lowlow <- divergent_txn_sort[which(divergent_txn_sort$left_gene_txn<120&divergent_txn_sort$right_gene_txn<120),]
nrow(highhigh);nrow(lowlow)

highhigh[,'width'] <- highhigh$right_coordinate-highhigh$left_coordinate+1
lowlow[,'width'] <- lowlow$right_coordinate-lowlow$left_coordinate+1
alldata <- list(lowlow$width,highhigh$width)
par(las=1)
boxplot(alldata,ylab="Intergenic size (bp)",names=c('low txn','high txn'),frame.plot=F,cex.lab=1.5,cex.axis=1.25,outline=F)
wilcox.test(lowlow$width,highhigh$width,paired=F) # p-value = 0.5389

# Pre-meiotic

txn_data <- 'X0hr.totRNA.rpkm'
transcription <- mrna0h
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
highhigh <- divergent_txn_sort[which(divergent_txn_sort$left_gene_txn>=40&divergent_txn_sort$right_gene_txn>=40),]
lowlow <- divergent_txn_sort[which(divergent_txn_sort$left_gene_txn<5&divergent_txn_sort$right_gene_txn<5),]
nrow(highhigh);nrow(lowlow)

highhigh[,'width'] <- highhigh$right_coordinate-highhigh$left_coordinate+1
lowlow[,'width'] <- lowlow$right_coordinate-lowlow$left_coordinate+1
alldata <- list(lowlow$width,highhigh$width)
par(las=1)
boxplot(alldata,ylab="Intergenic size (bp)",names=c('low txn','high txn'),frame.plot=F,cex.lab=1.5,cex.axis=1.25,outline=F)
wilcox.test(lowlow$width,highhigh$width,paired=F) # p-value = 6.941e-08
