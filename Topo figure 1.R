# Topo figure 1
####################################################################################
####################################################################################
library(hwglabr2)
library(GenomicRanges)
library(EnrichedHeatmap)
library(circlize)
library(ggplot2)

# Figure 1a

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
  bins_samp1 <- GenomeInfoDb::keepSeqlevels(bins_samp1, paste0("chr",chrnum))
  positions_samp1 <- bins_samp1@ranges@start + floor(bins_samp1@ranges@width / 2)
  df_samp1 <- data.frame(position=positions_samp1/1000, signal=bins_samp1$binned_score)
  return(df_samp1)
}
# function to plot genome data
plot_genomeview <- function(df_samp1,df_samp2,position1,position2,name1,name2,chrnum,color1,color2) {
  par(las=1)
  par(mfrow=c(2,1))
  ax=c(min(df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,1]),df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,1],max(df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,1]))
  ay=c(0,df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,2],0)
  bx=c(min(df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,1]),df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,1],max(df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,1]))
  by=c(0,df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,2],0)
  plot(df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,],
       xlab=paste0('Position on chromosome ',chrnum,' (kb)'), ylab=name1, type='h',col=color1,frame.plot=F)#,
  polygon(ax,ay,col=color1,border = NA)
  plot(df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,],
       xlab=paste0('Position on chromosome ',chrnum,' (kb)'), ylab=name2, type='h',col=color2,frame.plot=F)#,
  polygon(bx,by,col=color2,border = NA)
  par(mfrow=c(1,1))

}

Top2_wtdgv <- genomeView(Top2_wtd,"II",10)
Top1_mycdgv <- genomeView(Top1_mycd,"II",10)

plot_genomeview(Top1_mycdgv,Top2_wtdgv,640,720,"Top1","Top2","II","green","blue")

####################################################################################
####################################################################################
# Figure 1b

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

gff <- hwglabr2::get_gff(genome = 'SK1Yue')

Top2signal_at_ORFs <- hwglabr2::signal_at_orf2(signal_data=Top2_wtd, gff=gff,
                                               write_to_file=FALSE)
Top2signal_at_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=Top2signal_at_ORFs,
                                                      ci=0.95, rep_bootstrap=1000,
                                                      na_rm=TRUE)
Top1signal_at_ORFs <- hwglabr2::signal_at_orf2(signal_data=Top1_mycd, gff=gff,
                                               write_to_file=FALSE)
Top1signal_at_metaORF <- hwglabr2::signal_mean_and_ci(signal_data=Top1signal_at_ORFs,
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
Top2_wt = import_bedGraph("Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
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

Top2_wt = import_bedGraph("Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
Top1_myc = import_bedGraph("AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_FE.bdg.gz")
mnase3 = rtracklayer::import.bedGraph("Nucleosome_reps-SK1-MACS2_treat_pileup.bdg")

gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}

Top2_wtd = gendiv(Top2_wt)
Top1_mycd = gendiv(Top1_myc)
mnase3d = gendiv(mnase3)

intergen <- hwglabr2::get_intergenic_regions("SK1Yue",as_gr=T)
prom <- intergen[intergen$type=="divergent"|intergen$type=='tandem']
mcols(prom)['widths'] <- width(prom)
prom <- prom[order(width(prom),decreasing = T)]
midpoint <- floor(width(prom) / 2)
start(prom) <- start(prom) + midpoint
end(prom) <- start(prom)
mcols(prom)['class'] <- DataFrame(class=c(rep(1:4, each=length(prom)/4),4,4))

# Fig 1d,f
prommat <- normalizeToMatrix(Top1_mycd, prom, value_column = "score",
                             extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
col_fun <- colorRamp2(quantile(prommat, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
partition <- prom$class
EnrichedHeatmap(prommat, col = col_fun, name = "Top1",
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:4),
                                                                         show_error = TRUE)),
                top_annotation_height = unit(5, "cm"), row_title_rot = 0,
                axis_name = c("-1 kb", "promoters", "1 kb"),
                split=prom$class,
                row_order = 1:length(prom))+
  Heatmap(partition, col = structure(1:4, names = as.character(1:4)), name = "",row_order = 1:length(prom),
          show_row_names = FALSE, width = unit(5, "mm"))

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
