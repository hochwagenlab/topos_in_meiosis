# Topo figure 4
####################################################################################
####################################################################################
library(EnrichedHeatmap)
library(circlize)
library(hwglabr2)
library(GenomicRanges)

spo11oligo <- rtracklayer::import.bedGraph("/Volumes/LabShare/Jonna/Spo11_oligo_mapping/SK1Yue/Spo11oligo_WT1_SRR-clip-MACS2_extsize37/Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg")

intergen <- hwglabr2::get_intergenic_regions("SK1Yue",as_gr=T)
prom <- intergen[intergen$type=="divergent"|intergen$type=='tandem']
mcols(prom)['widths'] <- width(prom)
prom <- prom[order(width(prom),decreasing = T)]
midpoint <- floor(width(prom) / 2)
start(prom) <- start(prom) + midpoint
end(prom) <- start(prom)
mcols(prom)['class'] <- DataFrame(class=c(rep(1:4, each=length(prom)/4),4,4))

prommat <- normalizeToMatrix(spo11oligo, prom, value_column = "score",
                             extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
col_fun <- colorRamp2(quantile(prommat, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
partition <- prom$class
EnrichedHeatmap(prommat, col = col_fun, name = "Spo11",
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:4),
                                                                         show_error = TRUE)),
                top_annotation_height = unit(5, "cm"), row_title_rot = 0,
                axis_name = c("-1 kb", "promoters", "1 kb"),
                split=prom$class,
                row_order = 1:length(prom))+
  Heatmap(partition, col = structure(1:4, names = as.character(1:4)), name = "",row_order = 1:length(prom),
          show_row_names = FALSE, width = unit(5, "mm"))
