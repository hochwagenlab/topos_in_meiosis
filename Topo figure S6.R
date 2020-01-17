# Topo figure S6
####################################################################################
####################################################################################
# Figure S6a

Top2_34 = hwglabr2::import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-34C-493-533-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-34C-493-533-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
top2_1 = hwglabr2::import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-top2-4-496-537-Reps-SK1Yue-B3W3-MACS2/Top2-top2-4-496-537-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")

gendiv = function(bdg) {
  gavg = hwglabr2::average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  print(mean(bdg_new$score))
  return(bdg_new)
}
Top2_top2nf = gendiv(top2_1)
Top2_34cnf = gendiv(Top2_34)

# get gene seq info
genome_info <- hwglabr2::get_chr_coordinates(genome='SK1Yue')

# Sort sequences and levels to make sure they match
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)

Top2_34cnf_signal <- sort(GenomeInfoDb::sortSeqlevels(Top2_34cnf))
Top2_top2nf_signal <- sort(GenomeInfoDb::sortSeqlevels(Top2_top2nf))

# Add info to signal object
GenomeInfoDb::seqlengths(Top2_34cnf_signal) <- GenomeInfoDb::seqlengths(genome_info)
GenomeInfoDb::seqlengths(Top2_top2nf_signal) <- GenomeInfoDb::seqlengths(genome_info)

bins1000 <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(Top2_34cnf_signal),
                                      tilewidth=1000,
                                      cut.last.tile.in.chrom=TRUE)

# Get signal as "RleList"; the signal is stored in the "score" metadata column
Top2_34cnf_signal_score <- GenomicRanges::coverage(Top2_34cnf_signal, weight="score")
Top2_top2nf_signal_score <- GenomicRanges::coverage(Top2_top2nf_signal, weight="score")

# Get signal per tile
Top2_34cnf_bins1000 <- GenomicRanges::binnedAverage(bins1000, Top2_34cnf_signal_score, "binned1000_score")
Top2_top2nf_bins1000 <- GenomicRanges::binnedAverage(bins1000, Top2_top2nf_signal_score, "binned1000_score")

Top2_34cnf_bins1000[Top2_34cnf_bins1000$binned1000_score == 0]$binned1000_score <- NA
Top2_top2nf_bins1000[Top2_top2nf_bins1000$binned1000_score == 0]$binned1000_score <- NA

datamatch <- findOverlaps(Top2_34cnf_bins1000,Top2_top2nf_bins1000)
datamatch <- data.frame(Top2_34cnf_bins1000[queryHits(datamatch)],Top2_top2nf_bins1000[subjectHits(datamatch)])

data <- data.frame(Top2_34cnf_bins1000$binned1000_score,Top2_top2nf_bins1000$binned1000_score)
data_rm <- data.frame(data[which(data$Top2_34cnf_bins1000.binned1000_score!='NA'),])
data_rm <- data.frame(data_rm[which(data_rm$Top2_top2nf_bins1000.binned1000_score!='NA'),])
plot(data_rm,ylim=c(0,4),xlim=c(0,4))
cor(x=Top2_34cnf_bins1000$binned1000_score, y=Top2_top2nf_bins1000$binned1000_score, use="complete.obs", method="pearson")
#[1] 0.6187444

####################################################################################
# Figure S6b
library(GenomicRanges)
library(ggplot2)

Top2_wt = hwglabr2::import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-wildtype-413-504-Reps-SK1Yue-B3W3-MACS2/Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
Top2_top2 = hwglabr2::import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1Yue_MACS2_FE/Top2-top2-4-496-537-Reps-SK1Yue-B3W3-MACS2/Top2-top2-4-496-537-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")

gendiv = function(bdg) {
  gavg = hwglabr2::average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top2_wtd = gendiv(Top2_wt)
Top2_top2d = gendiv(Top2_top2)

gff <- hwglabr2::get_gff(genome = 'SK1Yue')

Top2signal_at_ORFs <- hwglabr2::signal_at_orf2(signal_data=Top2_wtd, gff=gff,
                                               write_to_file=FALSE)
Top2signal_at_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=Top2signal_at_ORFs,
                                                      ci=0.95, rep_bootstrap=1000,
                                                      na_rm=TRUE)
top2_1d_at_ORFs <- hwglabr2::signal_at_orf2(signal_data=Top2_top2d, gff=gff,
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