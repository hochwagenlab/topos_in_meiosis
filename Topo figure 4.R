# Topo figure 4
####################################################################################
####################################################################################
library(EnrichedHeatmap)
library(circlize)
library(hwglabr2)
library(GenomicRanges)
# Fig 4a
## Spo11 oligos lined up at Start, sorted by txn level
spo11oligo <- rtracklayer::import.bedGraph("/Volumes/LabShare/Jonna/Spo11_oligo_mapping/SK1Yue/Spo11oligo_WT1_SRR-clip-MACS2_extsize37/Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg")
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
mcols(gffall_sort) <- DataFrame(class=c(rep(1:4, each=length(gffall_sort)/4),4,4))

spo11oligo_atg <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, gffall_sort,
                                                     extend=c(500,500), w=10,empty_value=NA,
                                                     mean_mode="weighted",
                                                     value_column='score')
col_fun <- colorRamp2(quantile(spo11oligo_atg, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
partition <- gffall_sort$class
EnrichedHeatmap(spo11oligo_atg, col = col_fun, name = "Spo11", row_title_rot = 0,
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:4),
                                                                         show_error = T)),
                top_annotation_height = unit(2, "cm"),
                row_order = 1:length(gffall_sort),
                split=gffall_sort$class,
                axis_name = c("-500", "ATG","500"))+
  Heatmap(partition, col = structure(1:4, names = as.character(1:4)), name = "",row_order = 1:length(gffall_sort),
          show_row_names = FALSE, width = unit(5, "mm"))


## Spo11 oligos lined up at STOPs, sorted by txn level
spo11oligo <- rtracklayer::import.bedGraph("/Volumes/LabShare/Jonna/Spo11_oligo_mapping/SK1Yue/Spo11oligo_WT1_SRR-clip-MACS2_extsize37/Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg")
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
negstrandswitch <- GRanges(seqnames = seqnames(negstrand),ranges = IRanges(start=start(negstrand),end=start(negstrand)),strand=strand(negstrand),transcription=negstrand$transcription,gene=negstrand$gene)
start(posstrand) <- end(posstrand)
gffall <- c(posstrand,negstrandswitch)
gffall_sort <- gffall[order(gffall$transcription,decreasing = T)]
gffall_sort <- gffall_sort[which(gffall_sort$transcription!='<NA>')]
mcols(gffall_sort) <- DataFrame(class=c(rep(1:4, each=length(gffall_sort)/4),4,4))

spo11oligo_stop <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, gffall_sort,
                                                      extend=c(500,500), w=10,empty_value=NA,
                                                      mean_mode="weighted",
                                                      value_column='score')
col_fun <- colorRamp2(quantile(spo11oligo_atg, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
partition <- gffall_sort$class
EnrichedHeatmap(spo11oligo_stop, col = col_fun, name = "DSBs", row_title_rot = 0,
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
# Fig 4b

spo11oligo <- rtracklayer::import.bedGraph("/Volumes/LabShare/Jonna/Spo11_oligo_mapping/SK1Yue/Spo11oligo_WT1_SRR-clip-MACS2_extsize37/Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg")

intergen <- hwglabr2::get_intergenic_regions("SK1Yue",as_gr=T)
prom <- intergen[intergen$type=="divergent"|intergen$type=='tandem']
prom <- prom[order(width(prom),decreasing = T)]
midpoint <- floor(width(prom) / 2)
start(prom) <- start(prom) + midpoint
end(prom) <- start(prom)

prommat <- normalizeToMatrix(spo11oligo, prom, value_column = "score",
                             extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
col_fun <- colorRamp2(quantile(prommat, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
EnrichedHeatmap(prommat, col = col_fun, name = "Spo11",
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = "purple4"),
                                                                         show_error = TRUE)),
                top_annotation_height = unit(2, "cm"), row_title_rot = 0,
                axis_name = c("-1 kb", "promoters", "1 kb"),
                row_order = 1:length(prom))
####################################################################################
####################################################################################
# Fig 4c

mnase3 = rtracklayer::import.bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Nucleosome_replicates_PM-MACS2/Nucleosome_reps-SK1-MACS2_treat_pileup.bdg")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
mnase3d = gendiv(mnase3)
intergen <- hwglabr2::get_intergenic_regions("SK1Yue",as_gr=T)
prom <- intergen[intergen$type=="divergent"|intergen$type=='tandem']
prom <- prom[order(width(prom),decreasing = T)]
midpoint <- floor(width(prom) / 2)
start(prom) <- start(prom) + midpoint
end(prom) <- start(prom)

prommat <- normalizeToMatrix(mnase3d, prom, value_column = "score",
                             extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
col_fun <- colorRamp2(quantile(prommat, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
EnrichedHeatmap(prommat, col = col_fun, name = "Nuc 3h",
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = "purple4"),
                                                                         show_error = TRUE)),
                top_annotation_height = unit(2, "cm"), row_title_rot = 0,
                axis_name = c("-1 kb", "promoters", "1 kb"),
                row_order = 1:length(prom))
####################################################################################
####################################################################################
# Figure 4d: 

# Top1 at hotspots
Top1_myc = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/AH9847Myc-3h-735-841-Reps-SK1Yue-B3W4-MACS2/AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top1_mycd = gendiv(Top1_myc)
SK1Yue_Spo11_DSBs <- get_dsb_hotspots('SK1Yue')
midpoint <- floor(GenomicRanges::width(SK1Yue_Spo11_DSBs) / 2)
GenomicRanges::start(SK1Yue_Spo11_DSBs) <- GenomicRanges::start(SK1Yue_Spo11_DSBs) + midpoint
GenomicRanges::end(SK1Yue_Spo11_DSBs) <- GenomicRanges::start(SK1Yue_Spo11_DSBs)
SK1Yue_Spo11_DSBs <- SK1Yue_Spo11_DSBs[order(SK1Yue_Spo11_DSBs$score), ]
hotspots_sorted <- sort(SK1Yue_Spo11_DSBs, by =~ score, decreasing = T)
mcols(hotspots_sorted) <- DataFrame(class=c(rep(1:5, each=length(hotspots_sorted)/5),5))

Top1d_at_hotspots <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, hotspots_sorted,
                                                        extend=1000, w=10,
                                                        mean_mode="weighted",
                                                        value_column="score")
col_fun <- colorRamp2(quantile(Top1d_at_hotspots, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
partition <- hotspots_sorted$class
EnrichedHeatmap(Top1d_at_hotspots, col = col_fun, name = "Top1", row_title_rot = 0,
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:5),
                                                                         show_error = TRUE)),
                top_annotation_height = unit(4, "cm"),
                row_order = 1:length(hotspots_sorted),
                split=hotspots_sorted$class,
                axis_name = c("-1 kb", "DSB hotspots", "1 kb"))+
  Heatmap(partition, col = structure(1:5, names = as.character(1:5)), name = "",row_order = 1:length(hotspots_sorted),
          show_row_names = FALSE, width = unit(2, "mm"))


# Top2 at hotspots
Top2_wt = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-413-504-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top2_wtd = gendiv(Top2_wt)
SK1Yue_Spo11_DSBs <- get_dsb_hotspots('SK1Yue')
midpoint <- floor(GenomicRanges::width(SK1Yue_Spo11_DSBs) / 2)
GenomicRanges::start(SK1Yue_Spo11_DSBs) <- GenomicRanges::start(SK1Yue_Spo11_DSBs) + midpoint
GenomicRanges::end(SK1Yue_Spo11_DSBs) <- GenomicRanges::start(SK1Yue_Spo11_DSBs)
SK1Yue_Spo11_DSBs <- SK1Yue_Spo11_DSBs[order(SK1Yue_Spo11_DSBs$score), ]
hotspots_sorted <- sort(SK1Yue_Spo11_DSBs, by =~ score, decreasing = T)
mcols(hotspots_sorted) <- DataFrame(class=c(rep(1:5, each=length(hotspots_sorted)/5),5))

Top2d_at_hotspots <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, hotspots_sorted,
                                                        extend=1000, w=10,
                                                        mean_mode="weighted",
                                                        value_column="score")
col_fun <- colorRamp2(quantile(Top2d_at_hotspots, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
partition <- hotspots_sorted$class
EnrichedHeatmap(Top2d_at_hotspots, col = col_fun, name = "Top2", row_title_rot = 0,
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:5),
                                                                         show_error = TRUE)),
                top_annotation_height = unit(4, "cm"),
                row_order = 1:length(hotspots_sorted),
                split=hotspots_sorted$class,
                axis_name = c("-1 kb", "DSB hotspots", "1 kb"))+
  Heatmap(partition, col = structure(1:5, names = as.character(1:5)), name = "",row_order = 1:length(hotspots_sorted),
          show_row_names = FALSE, width = unit(2, "mm"))

####################################################################################
####################################################################################
# Fig 4e
extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "numeric")
import.broadPeak <- function(...) {
  rtracklayer::import(..., format="BED", extraCols=extraCols_broadPeak)
}

SK1Yue_Spo11_DSB <- get_dsb_hotspots('SK1Yue')
Top2_peak = import.broadPeak("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-413-504-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_peaks.broadPeak")
Top1_peak = import.broadPeak("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/AH9847Myc-3h-735-841-Reps-SK1Yue-B3W4-MACS2/AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_peaks.broadPeak")

#Hotspots that overlap with Top2
Spo11Top2 = findOverlaps(Top2_peak,SK1Yue_Spo11_DSB)
Spo11Top2 = SK1Yue_Spo11_DSB[subjectHits(Spo11Top2)]
#Hotspots that overlap with Top1
Spo11Top1 = findOverlaps(Top1_peak,SK1Yue_Spo11_DSB)
Spo11Top1 = SK1Yue_Spo11_DSB[subjectHits(Spo11Top1)]

#Find Top2 hotspots that overlap with Top1
Spo11Top2Top1 = findOverlaps(Spo11Top1,Spo11Top2)
Spo11Top2Top1 = Spo11Top2[subjectHits(Spo11Top2Top1)] #n=1694

#find Top2 peaks that do not overlap with the rest
Top2alone = Spo11Top2[!(Spo11Top2 %over% Spo11Top1)] #n=780

#find Top1 peaks that do not overlap with the rest
Top1alone = Spo11Top1[!(Spo11Top1 %over% Spo11Top2)] #n=253

#hotspots that do not overlap with Top1 or Top2
SK1Yue_Spo11_DSBalone = SK1Yue_Spo11_DSB[!(SK1Yue_Spo11_DSB %over% Spo11Top1)]
SK1Yue_Spo11_DSBalone = SK1Yue_Spo11_DSBalone[!(SK1Yue_Spo11_DSBalone %over% Spo11Top2)] # n=464

alldata=list(SK1Yue_Spo11_DSBalone$score, #no overlaps
             Top1alone$score, #Top1 alone
             Top2alone$score, #Top2 alone
             Spo11Top2Top1$score) #Top1+Top2
boxplot(alldata,names=c("None", "Top1 Only", "Top2 Only","Top1&Top2"),outline=F,ylab="Signal Intensity",frame.plot=F,cex.lab=1.5,cex.axis=1.25)

p <- c(wilcox.test(SK1Yue_Spo11_DSBalone$score,Top1alone$score,paired = F)$p.value,
              wilcox.test(SK1Yue_Spo11_DSBalone$score,Top2alone$score,paired = F)$p.value,
              wilcox.test(SK1Yue_Spo11_DSBalone$score,Spo11Top2Top1$score,paired = F)$p.value,
              wilcox.test(Top1alone$score,Top2alone$score,paired = F)$p.value,
              wilcox.test(Top1alone$score,Spo11Top2Top1$score,paired = F)$p.value,
              wilcox.test(Top2alone$score,Spo11Top2Top1$score,paired = F)$p.value)
p.adjust(p, method = "bonferroni", n = length(p))
