# Topo figure S1
####################################################################################
####################################################################################
# Figure S1c
# Top1
library(EnrichedHeatmap)
library(hwglabr2)
#convergent regions
Top1_myc = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/AH9847Myc-3h-735-841-Reps-SK1Yue-B3W4-MACS2/AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top1_mycd = gendiv(Top1_myc)
intergen = hwglabr2::get_intergenic_regions("SK1Yue")
conv = intergen[intergen$type == "convergent",]
midpoint = floor((conv$right_coordinate + conv$left_coordinate) / 2)
convmidpt <- GRanges(seqnames = conv$chr,ranges = IRanges(midpoint,midpoint))

Top1d_conv <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, convmidpt,
                                                 extend=1000, w=1,
                                                 mean_mode="weighted",
                                                 value_column="score")
Top1d_convci <- hwglabr2::signal_mean_and_ci(signal_data=Top1d_conv,
                                             ci=0.95, rep_bootstrap=1000,
                                             na_rm=TRUE)

# divergent regions
intergen = hwglabr2::get_intergenic_regions("SK1Yue")
div = intergen[intergen$type == "divergent",]
dmidpoint = floor((div$right_coordinate + div$left_coordinate) / 2)
divmidpt <- GRanges(seqnames = div$chr,ranges = IRanges(dmidpoint,dmidpoint))

Top1d_div <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, divmidpt,
                                                extend=1000, w=1,
                                                mean_mode="weighted",
                                                value_column="score")
Top1d_divci <- hwglabr2::signal_mean_and_ci(signal_data=Top1d_div,
                                            ci=0.95, rep_bootstrap=1000,
                                            na_rm=TRUE)

# tandem regions
intergen = hwglabr2::get_intergenic_regions("SK1Yue")
tand = intergen[intergen$type == "tandem",]
tmidpoint = floor((tand$right_coordinate + tand$left_coordinate) / 2)
tandmidpt <- GRanges(seqnames = tand$chr,ranges = IRanges(tmidpoint,tmidpoint))

Top1d_tand <- EnrichedHeatmap::normalizeToMatrix(Top1_mycd, tandmidpt,
                                                 extend=1000, w=1,
                                                 mean_mode="weighted",
                                                 value_column="score")
Top1d_tandci <- hwglabr2::signal_mean_and_ci(signal_data=Top1d_tand,
                                             ci=0.95, rep_bootstrap=1000,
                                             na_rm=TRUE)

group1s_gg <- data.frame(Data="convergent",Position=seq(1, 2000), Top1d_convci)
group2s_gg <- data.frame(Data="divergent",Position=seq(1, 2000), Top1d_divci)
group3s_gg <- data.frame(Data="tandem",Position=seq(1, 2000), Top1d_tandci)
allgroups <- rbind(group1s_gg,group2s_gg,group3s_gg)
# Set up the plot
library(ggplot2)
p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = 1000, lty = 3) +
  scale_x_continuous(breaks = c(0, 1000, 2000),
                     labels = c('-1kb', 'midpoint','1kb'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p