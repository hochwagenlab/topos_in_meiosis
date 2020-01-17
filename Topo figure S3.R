# Topo figure S3
####################################################################################
####################################################################################
library(hwglabr2)
library(GenomicRanges)
library(EnrichedHeatmap)
library(circlize)
library(ggplot2)


# Figure S3a
Top2_0 = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-0h-412-503-530-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-0h-412-503-530-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
nf_Top2_0h = 0.7541135 # calculated in Fig 3
Top2_3 = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-413-504-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
nf_Top2_3h = 1.013089 # calculated in Fig 3
median_nf = function(bdg,nf) {
  gmedian = median(rep(GenomicRanges::score(bdg),GenomicRanges::width(bdg)))
  print(gmedian)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score-gmedian
  print(mean(bdg_new$score))
  bdg_new$score <- (bdg_new$score)*nf
  print(mean(bdg_new$score))
  return(bdg_new)
}
Top2_0d = median_nf(Top2_0,nf_Top2_0h)
Top2_3d = median_nf(Top2_3,nf_Top2_3h)

intergen <- get_intergenic_regions('SK1Yue',as_gr=T)
promoter <- intergen[intergen$type=="divergent"|intergen$type=="tandem"]
mcols(promoter)['widths'] <- width(promoter)
promoter <- promoter[order(promoter$widths,decreasing=T)]
mcols(promoter)['class'] <- DataFrame(class=c(rep(1:3, each=length(promoter)/3),3,3))
midpoint <- floor(width(promoter) / 2)
start(promoter) <- start(promoter) + midpoint
end(promoter) <- start(promoter)

# 0h data
prom1 <- normalizeToMatrix(Top2_0d, promoter[promoter$class==1], value_column = "score",
                           extend = 1000, mean_mode = "weighted", w = 1,empty_value=NA)
prom2 <- normalizeToMatrix(Top2_0d, promoter[promoter$class==2], value_column = "score",
                           extend = 1000, mean_mode = "weighted", w = 1,empty_value=NA)
prom3 <- normalizeToMatrix(Top2_0d, promoter[promoter$class==3], value_column = "score",
                           extend = 1000, mean_mode = "weighted", w = 1,empty_value=NA)
prom1ci <- hwglabr2::signal_mean_and_ci(signal_data=prom1,
                                        ci=0.95, rep_bootstrap=1000,
                                        na_rm=TRUE)
prom2ci <- hwglabr2::signal_mean_and_ci(signal_data=prom2,
                                        ci=0.95, rep_bootstrap=1000,
                                        na_rm=TRUE)
prom3ci <- hwglabr2::signal_mean_and_ci(signal_data=prom3,
                                        ci=0.95, rep_bootstrap=1000,
                                        na_rm=TRUE)
prom1df0 <- data.frame(Data='Wide',Position=seq(-999, 1000), prom1ci)
prom2df0 <- data.frame(Data='Mid',Position=seq(-999, 1000), prom2ci)
prom3df0 <- data.frame(Data='Narrow',Position=seq(-999, 1000), prom3ci)
allgroup0 <- rbind(prom1df0,prom2df0,prom3df0)

p <- ggplot(allgroup0, aes(Position, Mean, group=Data, fill=Data,color=Data)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "promoters", y = "Top2 0h\nChIP-seq signal") +
  geom_vline(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-999,0,1000),
                     labels = c("-1kb","midpoint","1kb")) +
  scale_y_continuous(breaks = c(0,0.3,0.6,0.9),labels = c(0,0.3,0.6,0.9),limits = c(0,0.9))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) + geom_line()
p

# 3 hour data
prom13 <- normalizeToMatrix(Top2_3d, promoter[promoter$class==1], value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 1,empty_value=NA)
prom23 <- normalizeToMatrix(Top2_3d, promoter[promoter$class==2], value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 1,empty_value=NA)
prom33 <- normalizeToMatrix(Top2_3d, promoter[promoter$class==3], value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 1,empty_value=NA)
prom1ci3 <- hwglabr2::signal_mean_and_ci(signal_data=prom13,
                                         ci=0.95, rep_bootstrap=1000,
                                         na_rm=TRUE)
prom2ci3 <- hwglabr2::signal_mean_and_ci(signal_data=prom23,
                                         ci=0.95, rep_bootstrap=1000,
                                         na_rm=TRUE)
prom3ci3 <- hwglabr2::signal_mean_and_ci(signal_data=prom33,
                                         ci=0.95, rep_bootstrap=1000,
                                         na_rm=TRUE)
prom1df3 <- data.frame(Data='Wide',Position=seq(-999, 1000), prom1ci3)
prom2df3 <- data.frame(Data='Mid',Position=seq(-999, 1000), prom2ci3)
prom3df3 <- data.frame(Data='Narrow',Position=seq(-999, 1000), prom3ci3)
allgroup3 <- rbind(prom1df3,prom2df3,prom3df3)

p <- ggplot(allgroup3, aes(Position, Mean, group=Data, fill=Data,color=Data)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "promoters", y = "Top2 3h\nChIP-seq signal") +
  geom_vline(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-999,0,1000),
                     labels = c("-1kb","midpoint","1kb")) +
  scale_y_continuous(breaks = c(0,0.3,0.6,0.9),labels = c(0,0.3,0.6,0.9),limits = c(0,0.9))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) + geom_line()
p

####################################################################################
####################################################################################
# Figure S3b

spo11oligo <- rtracklayer::import.bedGraph("/Volumes/LabShare/Jonna/Spo11_oligo_mapping/SK1Yue/Spo11oligo_WT1_SRR-clip-MACS2_extsize37/Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg")

gff <- hwglabr2::get_gff('SK1Yue')
transcription <- read.csv('/Volumes/LabShare/HTGenomics/HiSeqOutputs/RNA-seq/2016.03.16-2h+3h/2017.06.16_SK1Yue_EdgeR_tpm.csv')
gff <- gff[which(gff$type=='gene')]
colnames(transcription)[1] <- "ID"
gff <- data.frame(gff)
gff_txn <- merge(x=gff,y=transcription[,c(1,3)],by='ID', all.x = TRUE)
gff_txn <- gff_txn[which(gff_txn$seqnames!='chrMT'),]
gff_txn <- gff_txn[which(gff_txn$seqnames!='scplasm1'),]

# Spo11 oligo signal at promoters
gffgr <- GRanges(seqnames = gff_txn$seqnames,IRanges(start = gff_txn$start,end=gff_txn$end),strand = gff_txn$strand,transcription=gff_txn$AH119_3h,gene=gff_txn$Name)
negstrand=gffgr[strand(gffgr)=='-']
posstrand=gffgr[strand(gffgr)=='+']
negstrandswitch <- GRanges(seqnames = seqnames(negstrand),ranges = IRanges(start=end(negstrand),end=end(negstrand)),strand=strand(negstrand),transcription=negstrand$transcription,gene=negstrand$gene)
end(posstrand) <- start(posstrand)
gffall <- c(posstrand,negstrandswitch)
gffall_sort <- gffall[order(gffall$transcription,decreasing = T)]
gffall_sort <- gffall_sort[which(gffall_sort$transcription!='<NA>')]

spo11oligo_atg <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, gffall_sort,
                                                     extend=c(500,0), w=1,empty_value=NA,
                                                     mean_mode="weighted",
                                                     value_column='score')
par(las=1)
sigrna <- data.frame(cbind(gffall_sort$transcription,rowSums(data.frame(spo11oligo_atg))))
colnames(sigrna) <- c('txn','spo11oligo')
plot(log(sigrna$txn),log(sigrna$spo11oligo))
plotline <- sigrna[which(sigrna$spo11oligo!=0),]
abline(lm(data=plotline,log(spo11oligo)~log(txn)))
cor(x=sigrna$txn,y=sigrna$spo11oligo, use="complete.obs", method="pearson")
#[1] 0.04665703
summary(lm(data=plotline,log(spo11oligo)~log(txn)))
####################################################################################
####################################################################################
# Figure S3c

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
                row_order = 1:length(hotspots_sorted),
                split=hotspots_sorted$class,
                axis_name = c("-1 kb", "DSB hotspots", "1 kb"))+
  Heatmap(partition, col = structure(1:5, names = as.character(1:5)), name = "",row_order = 1:length(hotspots_sorted),
          show_row_names = FALSE, width = unit(2, "mm"))

# average lines
Top1d_1 <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, hotspots_sorted[which(hotspots_sorted$class == 1)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top1d_1ci <- hwglabr2::signal_mean_and_ci(signal_data=Top1d_1,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
Top1d_2 <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, hotspots_sorted[which(hotspots_sorted$class == 2)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top1d_2ci <- hwglabr2::signal_mean_and_ci(signal_data=Top1d_2,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
Top1d_3 <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, hotspots_sorted[which(hotspots_sorted$class == 3)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top1d_3ci <- hwglabr2::signal_mean_and_ci(signal_data=Top1d_3,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
Top1d_4 <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, hotspots_sorted[which(hotspots_sorted$class == 4)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top1d_4ci <- hwglabr2::signal_mean_and_ci(signal_data=Top1d_4,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
Top1d_5 <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, hotspots_sorted[which(hotspots_sorted$class == 5)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top1d_5ci <- hwglabr2::signal_mean_and_ci(signal_data=Top1d_5,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)

group1_gg <- data.frame(Data="1",Position=seq(1, 2000), Top1d_1ci)
group2_gg <- data.frame(Data="2",Position=seq(1, 2000), Top1d_2ci)
group3_gg <- data.frame(Data="3",Position=seq(1, 2000), Top1d_3ci)
group4_gg <- data.frame(Data="4",Position=seq(1, 2000), Top1d_4ci)
group5_gg <- data.frame(Data="5",Position=seq(1, 2000), Top1d_5ci)
allgroups <- rbind(group1_gg,group2_gg,group3_gg,group4_gg,group5_gg)

# Set up the plot
p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = 1000, lty = 3) +
  scale_x_continuous(breaks = c(0, 1000, 2000),
                     labels = c('-1kb', 'midpoint','1kb'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p
######
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
                row_order = 1:length(hotspots_sorted),
                split=hotspots_sorted$class,
                axis_name = c("-1 kb", "DSB hotspots", "1 kb"))+
  Heatmap(partition, col = structure(1:5, names = as.character(1:5)), name = "",row_order = 1:length(hotspots_sorted),
          show_row_names = FALSE, width = unit(2, "mm"))

# average lines
Top2d_1 <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, hotspots_sorted[which(hotspots_sorted$class == 1)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top2d_1ci <- hwglabr2::signal_mean_and_ci(signal_data=Top2d_1,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
Top2d_2 <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, hotspots_sorted[which(hotspots_sorted$class == 2)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top2d_2ci <- hwglabr2::signal_mean_and_ci(signal_data=Top2d_2,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
Top2d_3 <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, hotspots_sorted[which(hotspots_sorted$class == 3)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top2d_3ci <- hwglabr2::signal_mean_and_ci(signal_data=Top2d_3,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
Top2d_4 <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, hotspots_sorted[which(hotspots_sorted$class == 4)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top2d_4ci <- hwglabr2::signal_mean_and_ci(signal_data=Top2d_4,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
Top2d_5 <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, hotspots_sorted[which(hotspots_sorted$class == 5)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top2d_5ci <- hwglabr2::signal_mean_and_ci(signal_data=Top2d_5,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)

group1_gg <- data.frame(Data="1",Position=seq(1, 2000), Top2d_1ci)
group2_gg <- data.frame(Data="2",Position=seq(1, 2000), Top2d_2ci)
group3_gg <- data.frame(Data="3",Position=seq(1, 2000), Top2d_3ci)
group4_gg <- data.frame(Data="4",Position=seq(1, 2000), Top2d_4ci)
group5_gg <- data.frame(Data="5",Position=seq(1, 2000), Top2d_5ci)
allgroups <- rbind(group1_gg,group2_gg,group3_gg,group4_gg,group5_gg)

p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = 1000, lty = 3) +
  scale_x_continuous(breaks = c(0, 1000, 2000),
                     labels = c('-1kb', 'midpoint','1kb'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p
####################################################################################
####################################################################################
# Figure S3d
spo11oligo <- rtracklayer::import.bedGraph("/Volumes/LabShare/Jonna/Spo11_oligo_mapping/SK1Yue/Spo11oligo_WT1_SRR-clip-MACS2_extsize37/Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg")

# Average signal plots by size
intergen <- hwglabr2::get_intergenic_regions("SK1Yue",as_gr=T)
prom <- intergen[intergen$type=="divergent"|intergen$type=='tandem']
mcols(prom)['widths'] <- width(prom)
prom <- prom[order(width(prom),decreasing = T)]
mcols(prom)['class'] <- DataFrame(class=c(rep(1:4, each=length(prom)/4),4,4))
midpoint <- floor(width(prom) / 2)
start(prom) <- start(prom) + midpoint
end(prom) <- start(prom)

spo11oligo_1 <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, prom[which(prom$class == 1)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
spo11oligo_1ci <- hwglabr2::signal_mean_and_ci(signal_data=spo11oligo_1,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
spo11oligo_2 <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, prom[which(prom$class == 2)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
spo11oligo_2ci <- hwglabr2::signal_mean_and_ci(signal_data=spo11oligo_2,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
spo11oligo_3 <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, prom[which(prom$class == 3)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
spo11oligo_3ci <- hwglabr2::signal_mean_and_ci(signal_data=spo11oligo_3,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
spo11oligo_4 <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, prom[which(prom$class == 4)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
spo11oligo_4ci <- hwglabr2::signal_mean_and_ci(signal_data=spo11oligo_4,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
group1_gg <- data.frame(Data="1",Position=seq(1, 2000), spo11oligo_1ci)
group2_gg <- data.frame(Data="2",Position=seq(1, 2000), spo11oligo_2ci)
group3_gg <- data.frame(Data="3",Position=seq(1, 2000), spo11oligo_3ci)
group4_gg <- data.frame(Data="4",Position=seq(1, 2000), spo11oligo_4ci)
allgroups <- rbind(group1_gg,group2_gg,group3_gg,group4_gg)

p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = 1000, lty = 3) +
  scale_x_continuous(breaks = c(0, 1000, 2000),
                     labels = c('-1kb', 'midpoint','1kb'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p

# Area under the curve
prommat1 <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, prom[which(prom$class==1)], value_column = "score",
                              extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
prommat2 <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, prom[which(prom$class==2)], value_column = "score",
                              extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
prommat3 <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, prom[which(prom$class==3)], value_column = "score",
                              extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
prommat4 <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, prom[which(prom$class==4)], value_column = "score",
                              extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
mean(rowSums(data.frame(prommat1)),na.rm=T)
#[1] 1192.41
mean(rowSums(data.frame(prommat2)),na.rm=T)
#[1] 926.3065
mean(rowSums(data.frame(prommat3)),na.rm=T)
#[1] 770.9869
mean(rowSums(data.frame(prommat4)),na.rm=T)
#[1] 827.0021

auc1data <- c(mean(rowSums(data.frame(prommat1)),na.rm=T),
              mean(rowSums(data.frame(prommat2)),na.rm=T),
              mean(rowSums(data.frame(prommat3)),na.rm=T),
              mean(rowSums(data.frame(prommat4)),na.rm=T))
auc1sd <- c(sd(rowSums(data.frame(prommat1)),na.rm=T)/sqrt(length(rowSums(data.frame(prommat1)))),
            sd(rowSums(data.frame(prommat2)),na.rm=T)/sqrt(length(rowSums(data.frame(prommat2)))),
            sd(rowSums(data.frame(prommat3)),na.rm=T)/sqrt(length(rowSums(data.frame(prommat3)))),
            sd(rowSums(data.frame(prommat4)),na.rm=T)/sqrt(length(rowSums(data.frame(prommat4)))))
prom1groups <- c('1','2','3','4')
auc1df <- data.frame(prom1groups,auc1data,auc1sd)

ggplot(auc1df,aes(x=prom1groups,y=auc1data,fill=prom1groups,width=0.8)) +
  stat_summary(fun.y=mean, geom="bar", width=0.5, alpha=0.25, colour=NA) +
  geom_errorbar(aes(ymin=auc1data-auc1sd, ymax=auc1data+auc1sd), width=.2,
                position=position_dodge(.9))+
  theme_classic()+
  labs(title = '', x = '', y = 'Relative Spo11 oligo amount')

wilcox.test(rowSums(data.frame(prommat1)),rowSums(data.frame(prommat4)))  
wilcox.test(rowSums(data.frame(prommat1)),rowSums(data.frame(prommat2)))  
wilcox.test(rowSums(data.frame(prommat2)),rowSums(data.frame(prommat3)))  
wilcox.test(rowSums(data.frame(prommat3)),rowSums(data.frame(prommat4)))  

p <- c(wilcox.test(rowSums(data.frame(prommat1)),rowSums(data.frame(prommat4)))$p.value,
       wilcox.test(rowSums(data.frame(prommat1)),rowSums(data.frame(prommat3)))$p.value,
       wilcox.test(rowSums(data.frame(prommat1)),rowSums(data.frame(prommat2)))$p.value,
       wilcox.test(rowSums(data.frame(prommat2)),rowSums(data.frame(prommat4)))$p.value,
       wilcox.test(rowSums(data.frame(prommat2)),rowSums(data.frame(prommat3)))$p.value,
       wilcox.test(rowSums(data.frame(prommat3)),rowSums(data.frame(prommat4)))$p.value)
p.adjust(p, method = "bonferroni", n = length(p))
#[1] 5.505085e-18 6.330337e-27 2.614342e-07 3.581659e-04 1.007698e-09 1.449708e-01

