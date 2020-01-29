# Topo figure S1
####################################################################################
####################################################################################
library(hwglabr2)
library(GenomicRanges)

# Figure S1C: correlation between Top1 and Top2

Top2_wt = import_bedGraph("Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
Top1_myc = import_bedGraph("AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz")
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

################################################################
################################################################
Top1_myc = import_bedGraph("AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
Top1_mycd = gendiv(Top1_myc)

# convergent regions
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
p <- ggplot(allgroups, aes(x=Position, y=Mean, group=Data, fill=Data,colour=Data))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_vline(xintercept = 1000, lty = 3) +
  scale_x_continuous(breaks = c(0, 1000, 2000),
                     labels = c('-1kb', 'midpoint','1kb'))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA) + geom_line()
p
