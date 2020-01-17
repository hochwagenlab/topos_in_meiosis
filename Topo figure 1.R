# Topo figure 1
####################################################################################
####################################################################################
library(hwglabr2)
library(GenomicRanges)
library(EnrichedHeatmap)
library(circlize)
library(ggplot2)

# Figure 1a

Top1_myc = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/AH9847Myc-3h-735-841-Reps-SK1Yue-B3W4-MACS2/AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz")
Top2_wt = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-413-504-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top1_mycd = gendiv(Top1_myc)
Top2_wtd = gendiv(Top2_wt)

# function to extract sequence from data sets
genomeView <- function(samp1,chrnum,tile) {
  genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
  samp1 <- sort(GenomeInfoDb::sortSeqlevels(samp1))
  GenomeInfoDb::seqlengths(samp1) <- GenomeInfoDb::seqlengths(genome_info)
  bins_samp1 <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(samp1),
                                          tilewidth=tile,
                                          cut.last.tile.in.chrom=TRUE)
  score_samp1 <- GenomicRanges::coverage(samp1, weight="score")
  bins_samp1 <- GenomicRanges::binnedAverage(bins_samp1, score_samp1, "binned_score")
  bins_samp1 <- GenomeInfoDb::keepSeqlevels(bins_samp1, paste0("chr",chrnum), pruning.mode="coarse")
  positions_samp1 <- bins_samp1@ranges@start + floor(bins_samp1@ranges@width / 2)
  df_samp1 <- data.frame(position=positions_samp1/1000, signal=bins_samp1$binned_score)
  return(df_samp1)
}
# function to plot genome data
plot_genomeview <- function(df_samp1,df_samp2,position1,position2,name1,name2,chrnum,color1,color2) {
  par(las=1)
  par(mfrow=c(2,1))
  ax=c(min(df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,1]),
       df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,1],
       max(df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,1]))
  ay=c(0,df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,2],0)
  bx=c(min(df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,1]),
       df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,1],
       max(df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,1]))
  by=c(0,df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,2],0)
  plot(df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,],
       xlab=paste0('Position on chromosome ',chrnum,' (kb)'), ylab=name1, 
       type='h',col=color1,frame.plot=F)#,
  polygon(ax,ay,col=color1,border = NA)
  plot(df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,],
       xlab=paste0('Position on chromosome ',chrnum,' (kb)'), ylab=name2, 
       type='h',col=color2,frame.plot=F)#,
  polygon(bx,by,col=color2,border = NA)
  par(mfrow=c(1,1))
  
}

Top1_mycdgv <- genomeView(Top1_mycd,"II",10)
Top2_wtdgv <- genomeView(Top2_wtd,"II",10)

plot_genomeview(Top1_mycdgv,Top2_wtdgv,640,720,"Top1","Top2","II","green","blue")

####################################################################################
####################################################################################
# Figure 1b

Top1_myc = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/AH9847Myc-3h-735-841-Reps-SK1Yue-B3W4-MACS2/AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz")
Top2_wt = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-413-504-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")

gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}

Top1_mycd = gendiv(Top1_myc)
Top2_wtd = gendiv(Top2_wt)

gff <- hwglabr2::get_gff(genome = 'SK1Yue')

Top1signal_at_ORFs <- hwglabr2::signal_at_orf2(signal_data=Top1_mycd, gff=gff,
                                               write_to_file=FALSE)
Top1signal_at_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=Top1signal_at_ORFs,
                                                      ci=0.95, rep_bootstrap=1000,
                                                      na_rm=TRUE)
Top2signal_at_ORFs <- hwglabr2::signal_at_orf2(signal_data=Top2_wtd, gff=gff,
                                               write_to_file=FALSE)
Top2signal_at_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=Top2signal_at_ORFs,
                                                      ci=0.95, rep_bootstrap=1000,
                                                      na_rm=TRUE)

group1s_gg <- data.frame(Data="Top1",Position=seq(1, 1000), Top1signal_at_metaORF)
group2s_gg <- data.frame(Data="Top2",Position=seq(1, 1000), Top2signal_at_metaORF)
allgroups <- rbind(group1s_gg,group2s_gg)

# Set up the plot
p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = c(250,750), lty = 3) +
  scale_x_continuous(breaks = c(250,750),
                     labels = c('start','stop'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p
####################################################################################
####################################################################################
# Figure 1c
# Top2
#convergent regions
Top2_wt = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-413-504-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top2_wtd = gendiv(Top2_wt)

# convergent regions
intergen = hwglabr2::get_intergenic_regions("SK1Yue")
conv = intergen[intergen$type == "convergent",]
midpoint = floor((conv$right_coordinate + conv$left_coordinate) / 2)
convmidpt <- GRanges(seqnames = conv$chr,ranges = IRanges(midpoint,midpoint))

Top2d_conv <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, convmidpt,
                                                 extend=1000, w=1,
                                                 mean_mode="weighted",
                                                 value_column="score")
Top2d_convci <- hwglabr2::signal_mean_and_ci(signal_data=Top2d_conv,
                                             ci=0.95, rep_bootstrap=1000,
                                             na_rm=TRUE)

# divergent regions
intergen = hwglabr2::get_intergenic_regions("SK1Yue")
div = intergen[intergen$type == "divergent",]
dmidpoint = floor((div$right_coordinate + div$left_coordinate) / 2)
divmidpt <- GRanges(seqnames = div$chr,ranges = IRanges(dmidpoint,dmidpoint))

Top2d_div <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, divmidpt,
                                                extend=1000, w=1,
                                                mean_mode="weighted",
                                                value_column="score")
Top2d_divci <- hwglabr2::signal_mean_and_ci(signal_data=Top2d_div,
                                            ci=0.95, rep_bootstrap=1000,
                                            na_rm=TRUE)

# tandem regions
intergen = hwglabr2::get_intergenic_regions("SK1Yue")
tand = intergen[intergen$type == "tandem",]
tmidpoint = floor((tand$right_coordinate + tand$left_coordinate) / 2)
tandmidpt <- GRanges(seqnames = tand$chr,ranges = IRanges(tmidpoint,tmidpoint))

Top2d_tand <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, tandmidpt,
                                                 extend=1000, w=1,
                                                 mean_mode="weighted",
                                                 value_column="score")
Top2d_tandci <- hwglabr2::signal_mean_and_ci(signal_data=Top2d_tand,
                                             ci=0.95, rep_bootstrap=1000,
                                             na_rm=TRUE)

group1s_gg <- data.frame(Data="convergent",Position=seq(1, 2000), Top2d_convci)
group2s_gg <- data.frame(Data="divergent",Position=seq(1, 2000), Top2d_divci)
group3s_gg <- data.frame(Data="tandem",Position=seq(1, 2000), Top2d_tandci)
allgroups <- rbind(group1s_gg,group2s_gg,group3s_gg)

# Set up the plot
p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = 1000, lty = 3) +
  scale_x_continuous(breaks = c(0, 1000, 2000),
                     labels = c('-1kb', 'midpoint','1kb'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p
####################################################################################
####################################################################################
# Figure 1d-g

Top2_wt = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-413-504-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
Top1_myc = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/AH9847Myc-3h-735-841-Reps-SK1Yue-B3W4-MACS2/AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz")

gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}

Top2_wtd = gendiv(Top2_wt)
Top1_mycd = gendiv(Top1_myc)

intergen <- hwglabr2::get_intergenic_regions("SK1Yue",as_gr=T)
prom <- intergen[intergen$type=="divergent"|intergen$type=='tandem']
mcols(prom)['widths'] <- width(prom)
prom <- prom[order(width(prom),decreasing = T)]
#mcols(prom)['class'] <- DataFrame(class=c(rep(1:4, each=length(prom)/4),4,4))
midpoint <- floor(width(prom) / 2)
start(prom) <- start(prom) + midpoint
end(prom) <- start(prom)

# Fig 1d
prommat <- normalizeToMatrix(Top1_mycd, prom, value_column = "score",
                             extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
col_fun <- colorRamp2(quantile(prommat, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
EnrichedHeatmap(prommat, col = col_fun, name = "Top1",
                axis_name = c("-1 kb", "promoters", "1 kb"),
                row_order = 1:length(prom))

# Fig 1f
# auc
prom1mat1 <- normalizeToMatrix(Top1_mycd, prom[which(prom$class==1)], value_column = "score",
                             extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
prom1mat2 <- normalizeToMatrix(Top1_mycd, prom[which(prom$class==2)], value_column = "score",
                             extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
prom1mat3 <- normalizeToMatrix(Top1_mycd, prom[which(prom$class==3)], value_column = "score",
                             extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
prom1mat4 <- normalizeToMatrix(Top1_mycd, prom[which(prom$class==4)], value_column = "score",
                             extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
mean(rowSums(data.frame(prom1mat1)),na.rm=T)
#[1] 230.1425
mean(rowSums(data.frame(prom1mat2)),na.rm=T)
#[1] 211.857
mean(rowSums(data.frame(prom1mat3)),na.rm=T)
#[1] 205.3683
mean(rowSums(data.frame(prom1mat4)),na.rm=T)
#[1] 199.506

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
  labs(title = '', x = '', y = 'Relative Top1 amount')

wilcox.test(rowSums(data.frame(prom1mat1)),rowSums(data.frame(prom1mat4)))  
wilcox.test(rowSums(data.frame(prom1mat1)),rowSums(data.frame(prom1mat2)))  
wilcox.test(rowSums(data.frame(prom1mat2)),rowSums(data.frame(prom1mat3)))  
wilcox.test(rowSums(data.frame(prom1mat3)),rowSums(data.frame(prom1mat4)))  
p <- c(wilcox.test(rowSums(data.frame(prom1mat1)),rowSums(data.frame(prom1mat4)))$p.value,
       wilcox.test(rowSums(data.frame(prom1mat1)),rowSums(data.frame(prom1mat3)))$p.value,
       wilcox.test(rowSums(data.frame(prom1mat1)),rowSums(data.frame(prom1mat2)))$p.value,
       wilcox.test(rowSums(data.frame(prom1mat2)),rowSums(data.frame(prom1mat4)))$p.value,
       wilcox.test(rowSums(data.frame(prom1mat2)),rowSums(data.frame(prom1mat3)))$p.value,
       wilcox.test(rowSums(data.frame(prom1mat3)),rowSums(data.frame(prom1mat4)))$p.value)
p.adjust(p, method = "bonferroni", n = length(p))
#[1] 9.000612e-111  3.679108e-74  4.666501e-39  6.875570e-43  1.622966e-12  1.307821e-11

# Average signal plots by size
intergen <- hwglabr2::get_intergenic_regions("SK1Yue",as_gr=T)
prom <- intergen[intergen$type=="divergent"|intergen$type=='tandem']
mcols(prom)['widths'] <- width(prom)
prom <- prom[order(width(prom),decreasing = T)]
mcols(prom)['class'] <- DataFrame(class=c(rep(1:4, each=length(prom)/4),4,4))
midpoint <- floor(width(prom) / 2)
start(prom) <- start(prom) + midpoint
end(prom) <- start(prom)

Top1d_1 <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, prom[which(prom$class == 1)],
                                                 extend=1000, w=1,
                                                 mean_mode="weighted",
                                                 value_column="score")
Top1d_1ci <- hwglabr2::signal_mean_and_ci(signal_data=Top1d_1,
                                             ci=0.95, rep_bootstrap=1000,
                                             na_rm=TRUE)
Top1d_2 <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, prom[which(prom$class == 2)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top1d_2ci <- hwglabr2::signal_mean_and_ci(signal_data=Top1d_2,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
Top1d_3 <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, prom[which(prom$class == 3)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top1d_3ci <- hwglabr2::signal_mean_and_ci(signal_data=Top1d_3,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
Top1d_4 <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, prom[which(prom$class == 4)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top1d_4ci <- hwglabr2::signal_mean_and_ci(signal_data=Top1d_4,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
group1_gg <- data.frame(Data="1",Position=seq(1, 2000), Top1d_1ci)
group2_gg <- data.frame(Data="2",Position=seq(1, 2000), Top1d_2ci)
group3_gg <- data.frame(Data="3",Position=seq(1, 2000), Top1d_3ci)
group4_gg <- data.frame(Data="4",Position=seq(1, 2000), Top1d_4ci)
allgroups <- rbind(group1_gg,group2_gg,group3_gg,group4_gg)

# Set up the plot
p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = 1000, lty = 3) +
  scale_x_continuous(breaks = c(0, 1000, 2000),
                     labels = c('-1kb', 'midpoint','1kb'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p
########################################################################################
# Fig 1e,g
prommat <- normalizeToMatrix(Top2_wtd, prom, value_column = "score",
                             extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
col_fun <- colorRamp2(quantile(prommat, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))

EnrichedHeatmap(prommat, col = col_fun, name = "Top2",
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:4),
                                                                         show_error = TRUE)),
                top_annotation_height = unit(5, "cm"), row_title_rot = 0,
                axis_name = c("-1 kb", "promoters", "1 kb"),
                split=prom$class,
                row_order = 1:length(prom))+
  Heatmap(partition, col = structure(1:4, names = as.character(1:4)), name = "",row_order = 1:length(prom),
          show_row_names = FALSE, width = unit(5, "mm"))

#auc
prom2mat1 <- normalizeToMatrix(Top2_wtd, prom[which(prom$class==1)], value_column = "score",
                              extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
prom2mat2 <- normalizeToMatrix(Top2_wtd, prom[which(prom$class==2)], value_column = "score",
                              extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
prom2mat3 <- normalizeToMatrix(Top2_wtd, prom[which(prom$class==3)], value_column = "score",
                              extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
prom2mat4 <- normalizeToMatrix(Top2_wtd, prom[which(prom$class==4)], value_column = "score",
                              extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
mean(rowSums(data.frame(prom2mat1)),na.rm=T)
#[1] 264.2913
mean(rowSums(data.frame(prom2mat2)),na.rm=T)
#[1] 227.5315
mean(rowSums(data.frame(prom2mat3)),na.rm=T)
#[1] 212.473
mean(rowSums(data.frame(prom2mat4)),na.rm=T)
#[1] 206.3432
aucdata <- c(mean(rowSums(data.frame(prom2mat1)),na.rm=T),
                        mean(rowSums(data.frame(prom2mat2)),na.rm=T),
                        mean(rowSums(data.frame(prom2mat3)),na.rm=T),
                        mean(rowSums(data.frame(prom2mat4)),na.rm=T))
aucsd <- c(sd(rowSums(data.frame(prom2mat1)),na.rm=T)/sqrt(length(rowSums(data.frame(prom2mat1)))),
           sd(rowSums(data.frame(prom2mat2)),na.rm=T)/sqrt(length(rowSums(data.frame(prom2mat2)))),
           sd(rowSums(data.frame(prom2mat3)),na.rm=T)/sqrt(length(rowSums(data.frame(prom2mat3)))),
           sd(rowSums(data.frame(prom2mat4)),na.rm=T)/sqrt(length(rowSums(data.frame(prom2mat4)))))
promgroups <- c('1','2','3','4')
aucdf <- data.frame(promgroups,aucdata,aucsd)

ggplot(aucdf,aes(x=promgroups,y=aucdata,fill=promgroups,width=0.8)) +
  stat_summary(fun.y=mean, geom="bar", width=0.5, alpha=0.25, colour=NA) +
  #geom_point(size=1.5, alpha=1) +
  geom_errorbar(aes(ymin=aucdata-aucsd, ymax=aucdata+aucsd), width=.2,
                position=position_dodge(.9))+
  theme_classic()+
  labs(title = '', x = '', y = 'Relative Top2 amount')

wilcox.test(rowSums(data.frame(prom2mat1)),rowSums(data.frame(prom2mat4)))  #p-value < 2.2e-16
wilcox.test(rowSums(data.frame(prom2mat1)),rowSums(data.frame(prom2mat2)))  #p-value = 4.357e-08
wilcox.test(rowSums(data.frame(prom2mat2)),rowSums(data.frame(prom2mat3)))  #p-value = 1.679e-10
wilcox.test(rowSums(data.frame(prom2mat3)),rowSums(data.frame(prom2mat4)))  #p-value = 0.02416
p <- c(wilcox.test(rowSums(data.frame(prom2mat1)),rowSums(data.frame(prom2mat4)))$p.value,
       wilcox.test(rowSums(data.frame(prom2mat1)),rowSums(data.frame(prom2mat3)))$p.value,
       wilcox.test(rowSums(data.frame(prom2mat1)),rowSums(data.frame(prom2mat2)))$p.value,
       wilcox.test(rowSums(data.frame(prom2mat2)),rowSums(data.frame(prom2mat4)))$p.value,
       wilcox.test(rowSums(data.frame(prom2mat2)),rowSums(data.frame(prom2mat3)))$p.value,
       wilcox.test(rowSums(data.frame(prom2mat3)),rowSums(data.frame(prom2mat4)))$p.value)
p.adjust(p, method = "bonferroni", n = length(p))
#[1] 5.571130e-75 2.847767e-55 7.776200e-24 3.165123e-28 2.461280e-12 8.520812e-05

# Average signal plots by size
intergen <- hwglabr2::get_intergenic_regions("SK1Yue",as_gr=T)
prom <- intergen[intergen$type=="divergent"|intergen$type=='tandem']
mcols(prom)['widths'] <- width(prom)
prom <- prom[order(width(prom),decreasing = T)]
mcols(prom)['class'] <- DataFrame(class=c(rep(1:4, each=length(prom)/4),4,4))
midpoint <- floor(width(prom) / 2)
start(prom) <- start(prom) + midpoint
end(prom) <- start(prom)

Top2d_1 <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, prom[which(prom$class == 1)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top2d_1ci <- hwglabr2::signal_mean_and_ci(signal_data=Top2d_1,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
Top2d_2 <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, prom[which(prom$class == 2)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top2d_2ci <- hwglabr2::signal_mean_and_ci(signal_data=Top2d_2,
                                           ci=0.95, rep_bootstrap=1000,
                                           na_rm=TRUE)
Top2d_3 <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, prom[which(prom$class == 3)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top2d_3ci <- hwglabr2::signal_mean_and_ci(signal_data=Top2d_3,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
Top2d_4 <- EnrichedHeatmap::normalizeToMatrix(Top2_wtd, prom[which(prom$class == 4)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
Top2d_4ci <- hwglabr2::signal_mean_and_ci(signal_data=Top2d_4,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
group1_gg <- data.frame(Data="1",Position=seq(1, 2000), Top2d_1ci)
group2_gg <- data.frame(Data="2",Position=seq(1, 2000), Top2d_2ci)
group3_gg <- data.frame(Data="3",Position=seq(1, 2000), Top2d_3ci)
group4_gg <- data.frame(Data="4",Position=seq(1, 2000), Top2d_4ci)
allgroups <- rbind(group1_gg,group2_gg,group3_gg,group4_gg)

# Set up the plot
p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = 1000, lty = 3) +
  scale_x_continuous(breaks = c(0, 1000, 2000),
                     labels = c('-1kb', 'midpoint','1kb'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p
########################################################################################
spo11oligo <- rtracklayer::import.bedGraph("/Volumes/LabShare/Jonna/Spo11_oligo_mapping/SK1Yue/Spo11oligo_WT1_SRR-clip-MACS2_extsize37/Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg")
intergen <- hwglabr2::get_intergenic_regions("SK1Yue",as_gr=T)
prom <- intergen[intergen$type=="divergent"|intergen$type=='tandem']
mcols(prom)['widths'] <- width(prom)
prom <- prom[order(width(prom),decreasing = T)]
mcols(prom)['class'] <- DataFrame(class=c(rep(1:4, each=length(prom)/4),4,4))
midpoint <- floor(width(prom) / 2)
start(prom) <- start(prom) + midpoint
end(prom) <- start(prom)

#auc
prommat1 <- normalizeToMatrix(spo11oligo, prom[which(prom$class==1)], value_column = "score",
                              extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
prommat2 <- normalizeToMatrix(spo11oligo, prom[which(prom$class==2)], value_column = "score",
                              extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
prommat3 <- normalizeToMatrix(spo11oligo, prom[which(prom$class==3)], value_column = "score",
                              extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
prommat4 <- normalizeToMatrix(spo11oligo, prom[which(prom$class==4)], value_column = "score",
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
  #geom_point(size=1.5, alpha=1) +
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
###################################
# correlation between Top1 and Top2

Top2_wt = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-413-504-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
Top1_myc = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/AH9847Myc-3h-735-841-Reps-SK1Yue-B3W4-MACS2/AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top2_wtd = gendiv(Top2_wt)
Top1_mycd = gendiv(Top1_myc)


# get gene seq info
genome_info <- hwglabr2::get_chr_coordinates(genome='SK1Yue')

# Sort sequences and levels to make sure they match
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)

Top2_wt_signal <- sort(GenomeInfoDb::sortSeqlevels(Top2_wt))
Top1_myc_signal <- sort(GenomeInfoDb::sortSeqlevels(Top1_myc))

# Add info to signal object
GenomeInfoDb::seqlengths(Top2_wt_signal) <- GenomeInfoDb::seqlengths(genome_info)
GenomeInfoDb::seqlengths(Top1_myc_signal) <- GenomeInfoDb::seqlengths(genome_info)

bins1000 <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(Top2_wt_signal),
                                  tilewidth=1000,
                                  cut.last.tile.in.chrom=TRUE)

# Get signal as "RleList"; the signal is stored in the "score" metadata column
Top2_wt_score <- GenomicRanges::coverage(Top2_wt_signal, weight="score")
Top1_score <- GenomicRanges::coverage(Top1_myc_signal, weight="score")

# Get signal per tile
Top2_wt_bins1000 <- GenomicRanges::binnedAverage(bins1000, Top2_wt_score, "binned1000_score")
Top1_bins1000 <- GenomicRanges::binnedAverage(bins1000, Top1_score, "binned1000_score")

Top2_wt_bins1000[Top2_wt_bins1000$binned1000_score == 0]$binned1000_score <- NA
Top1_bins1000[Top1_bins1000$binned1000_score == 0]$binned1000_score <- NA

datamatch <- findOverlaps(Top2_wt_bins1000,Top1_bins1000)
datamatch <- data.frame(Top2_wt_bins1000[queryHits(datamatch)],Top1_bins1000[subjectHits(datamatch)])

data <- data.frame(Top2_wt_bins1000$binned1000_score,Top1_bins1000$binned1000_score)
data_rm <- data.frame(data[which(data$Top2_wt_bins1000.binned1000_score!='NA'),])
data_rm <- data.frame(data_rm[which(data_rm$Top1_bins1000.binned1000_score!='NA'),])
plot(data_rm,ylim=c(0,3.25),xlim=c(0,3.25))
abline(lm(data=data_rm,Top1_bins1000.binned1000_score~Top2_wt_bins1000.binned1000_score))
cor(x=Top2_wt_bins1000$binned1000_score, y=Top1_bins1000$binned1000_score, use="complete.obs", method="pearson")
#[1] 0.7557327


