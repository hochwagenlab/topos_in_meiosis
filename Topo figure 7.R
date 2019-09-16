# Topo figure 6
####################################################################################
####################################################################################
library(hwglabr2)
library(GenomicRanges)
library(EnrichedHeatmap)
library(ggplot2)

# Figure 6B

spo11oligo <- rtracklayer::import.bedGraph("/Volumes/LabShare/Jonna/Spo11_oligo_mapping/SK1Yue/Spo11oligo_WT1_SRR-clip-MACS2_extsize37/Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg")
Top2_wt34 = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-34C-493-533-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-34C-493-533-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
Top2_top2 = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-top2-4-496-537-Reps-SK1Yue-B3W3-MACS2/Top2-top2-4-496-537-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
Red1_WT = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Red1-WT-34C-410-495-528-Reps-SK1Yue-B3W3-MACS2/Red1-WT-34C-410-495-528-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
Red1_top2 = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Red1-top2-4-411-498-535-Reps-SK1Yue-B3W3-MACS2/Red1-top2-4-411-498-535-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")

gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}

Red1_WTd = gendiv(Red1_WT)
Top2_wt34d = gendiv(Top2_wt34)
Top2_top2d = gendiv(Top2_top2)
Red1_top2d = gendiv(Red1_top2)

SK1Yue_Red1_summits <- get_Red1_summits("SK1Yue")
SK1Yue_Red1_summits <- SK1Yue_Red1_summits[order(SK1Yue_Red1_summits$score),]

genomeView <- function(samp1,chrnum,tile) {
  genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
  samp1 <- sort(GenomeInfoDb::sortSeqlevels(samp1))
  GenomeInfoDb::seqlengths(samp1) <- GenomeInfoDb::seqlengths(genome_info)
  bins_samp1 <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(samp1),
                                          tilewidth=tile,
                                          cut.last.tile.in.chrom=TRUE)
  score_samp1 <- GenomicRanges::coverage(samp1, weight="score")
  bins_samp1 <- GenomicRanges::binnedAverage(bins_samp1, score_samp1, "binned_score")
  bins_samp1 <- GenomeInfoDb::keepSeqlevels(bins_samp1, paste0("chr",chrnum))
  positions_samp1 <- bins_samp1@ranges@start + floor(bins_samp1@ranges@width / 2)
  df_samp1 <- data.frame(position=positions_samp1/1000, signal=bins_samp1$binned_score)
  return(df_samp1)
}

plot_genomeview <- function(df_samp1,df_samp2,position1,position2,name1,name2,chrnum,color1,color2) {
  par(las=1)
  par(mfrow=c(2,1))
  ax=c(min(df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,1]),df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,1],max(df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,1]))
  ay=c(0,df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,2],0)
  bx=c(min(df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,1]),df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,1],max(df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,1]))
  by=c(0,df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,2],0)
  plot(df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,], xlab=paste0('Position on chromosome ',chrnum,' (kb)'), ylab=name1, type='h',col=color1,frame.plot=F)
  polygon(ax,ay,col=color1,border = NA)
  plot(df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,], xlab=paste0('Position on chromosome ',chrnum,' (kb)'), ylab=name2, type='h',col=color2,frame.plot=F)
  polygon(bx,by,col=color2,border = NA)
  par(mfrow=c(1,1))
  
}

Top2_wtdgv <- genomeView(Top2_wt34d,"XII",10)
Top2_1gv <- genomeView(Top2_top2d,"XII",10)
Red1_WTgv <- genomeView(Red1_WTd,"XII",10)
spo11oligogv <- genomeView(spo11oligo,"XII",10)
Red1_top2gv <- genomeView(Red1_top2d,"XII",10)

plot_genomeview(Top2_wtdgv,Top2_1gv,103,155,"Top2","Top2-1","XII","blue","lightgreen")
plot_genomeview(spo11oligogv,Red1_WTgv,103,155,"Spo11","Red1","XII","black","red")
plot_genomeview(Red1_WTgv,Red1_top2gv,103,155,"Spo11","Red1","XII","red","red4")

###################################################
# Figure 6C

Top2_wt = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-34C-493-533-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-34C-493-533-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
top2_1 = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-top2-4-496-537-Reps-SK1Yue-B3W3-MACS2/Top2-top2-4-496-537-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")

gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}

Top2_wtd = gendiv(Top2_wt)
top2_1d = gendiv(top2_1)

gff <- hwglabr2::get_gff(genome = 'SK1Yue')

Top2signal_at_ORFs <- hwglabr2::signal_at_orf2(signal_data=Top2_wtd, gff=gff,
                                               write_to_file=FALSE)
Top2signal_at_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=Top2signal_at_ORFs,
                                                      ci=0.95, rep_bootstrap=1000,
                                                      na_rm=TRUE)
top2_1d_at_ORFs <- hwglabr2::signal_at_orf2(signal_data=top2_1d, gff=gff,
                                               write_to_file=FALSE)
top2_1d_at_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=top2_1d_at_ORFs,
                                                      ci=0.95, rep_bootstrap=1000,
                                                      na_rm=TRUE)

group1s_gg <- data.frame(Data="Top2",Position=seq(1, 1000), Top2signal_at_metaORF)
group2s_gg <- data.frame(Data="Top2-1",Position=seq(1, 1000), top2_1d_at_metaORF)
allgroups <- rbind(group1s_gg,group2s_gg)
# Set up the plot
p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = c(250,750), lty = 3) +
  scale_x_continuous(breaks = c(250,750),
                     labels = c('start','stop'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p

###################################################

# Figure 6D
#convergent regions
Top2_top2 = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-top2-4-496-537-Reps-SK1Yue-B3W3-MACS2/Top2-top2-4-496-537-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top2_top2d = gendiv(Top2_top2)

# convergent regions
intergen = hwglabr2::get_intergenic_regions("SK1Yue")
conv = intergen[intergen$type == "convergent",]
midpoint = floor((conv$right_coordinate + conv$left_coordinate) / 2)
convmidpt <- GRanges(seqnames = conv$chr,ranges = IRanges(midpoint,midpoint))

Top2d_conv <- EnrichedHeatmap::normalizeToMatrix(Top2_top2d, convmidpt,
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

Top2d_div <- EnrichedHeatmap::normalizeToMatrix(Top2_top2d, divmidpt,
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

Top2d_tand <- EnrichedHeatmap::normalizeToMatrix(Top2_top2d, tandmidpt,
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
                     labels = c('-1kb', 'midpoint','1kb')) +
  scale_y_continuous(breaks = c(0.6, 0.9,1.2,1.5,1.8),
                     limits=c(0.8,1.8))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p

