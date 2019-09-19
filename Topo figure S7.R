# Topo figure S6
####################################################################################
####################################################################################
library(hwglabr2)
library(GenomicRanges)
library(EnrichedHeatmap)
library(ggplot2)

# Figure 6c

Top2_wt = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-34C-493-533-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-34C-493-533-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
top2_1 = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-top2-4-496-537-Reps-SK1Yue-B3W3-MACS2/Top2-top2-4-496-537-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
Top2_HA = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/2014-10-13/SK1_Yue/AH7688B-141013-SK1Yue-B3W3-MACS2/AH7688B-141013-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")

gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}

Top2_wtd = gendiv(Top2_wt)
top2_1d = gendiv(top2_1)
Top2_HAd = gendiv(Top2_HA)

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
Top2HAsignal_at_ORFs <- hwglabr2::signal_at_orf2(signal_data=Top2_HAd, gff=gff,
                                                 write_to_file=FALSE)
Top2HAsignal_at_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=Top2HAsignal_at_ORFs,
                                                        ci=0.95, rep_bootstrap=1000,
                                                        na_rm=TRUE)

group1s_gg <- data.frame(Data="Top2",Position=seq(1, 1000), Top2signal_at_metaORF)
group2s_gg <- data.frame(Data="Top2-1",Position=seq(1, 1000), top2_1d_at_metaORF)
group3s_gg <- data.frame(Data="Top2-HA",Position=seq(1, 1000), Top2HAsignal_at_metaORF)
allgroups <- rbind(group1s_gg,group2s_gg,group3s_gg)
# Set up the plot
p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = c(250,750), lty = 3) +
  scale_x_continuous(breaks = c(250,750),
                     labels = c('start','stop'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p