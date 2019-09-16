# Topo figure 3
####################################################################################
####################################################################################
library(stringr)
library(readr)
library(hwglabr2)
library(ggplot2)
library(EnrichedHeatmap)

# Figure 3a

spikein_normalization_factor_from_counts <- function(
  ref_chip_counts, ref_input_counts, test_chip_counts, test_input_counts,
  return_counts=FALSE) {
  
  # Put paths in list
  files <- list(ref_chip=ref_chip_counts, ref_input=ref_input_counts,
                test_chip=test_chip_counts, test_input=test_input_counts)
  
  # Convert each element into list, if not one already
  for (i in seq_along(files)) {
    if (!is.list(files[[i]])) files[[i]] <- list(files[[i]])
  }
  
  # Print files to read to console
  message('>>> Read alignment count files:')
  for (i in seq_along(files)) {
    for (file in files[[i]]) {
      message('   ', basename(file))
    }
  }    
  
  message()
  # Read files into tibble in list
  tables <- list()
  for (i in seq_along(files)) {
    tables[[i]] <- sapply(files[[i]], FUN=read_tsv, col_names=F,
                          simplify=FALSE, USE.NAMES=TRUE)
  }
  names(tables) <- names(files)
  
  message()
  # Get read counts per chromosome
  message('>>> Count reads per genome:')
  counts <- list()
  for (i in seq_along(tables)) {
    counts[[i]] <- sapply(tables[[i]], FUN=sum_per_genome,
                          simplify=FALSE, USE.NAMES=TRUE)
  }
  names(counts) <- names(tables)
  
  # Add-up counts for replicates (results in nested lists)
  for (i in seq_along(counts)) {
    if (length(counts[[i]]) > 1) {
      total <- counts[[i]][[1]]
      for (j in 2:length(counts[[i]])) {
        total <- total + counts[[i]][[j]]
      }
      counts[[i]] <- total
    } else counts[[i]] <- unlist(counts[[i]])
  }
  
  if (return_counts) {
    message('---')
    message('Done!')
    return(counts)
  }
  
  # Compute normalization factor
  result <- normalization_factor(ctrl_input=counts$ref_input,
                                 ctrl_chip=counts$ref_chip,
                                 test_input=counts$test_input,
                                 test_chip=counts$test_chip)
  
  message('---')
  message('Done!')
  
  return(result)
}

# Helper functions
sum_per_genome <- function(df) {
  # Compute sum of reads aligned to each genome
  S288C <- sum(
    df[apply(df, 1, function(x) str_detect(x[1],'_S288C')), 2])
  SK1 <- sum(
    df[apply(df, 1, function(x) str_detect(x[1], '_SK1')), 2])
  
  # Print result to console
  message('  S288C: ', formatC(S288C, big.mark=",",
                               drop0trailing=TRUE, format="f"))
  message('  SK1: ', formatC(SK1, big.mark=",",
                             drop0trailing=TRUE, format="f"))
  message('      ', round(S288C * 100 / (SK1 + S288C), 1), '% spike-in reads')
  
  # Return result as named vector
  c('S288C'=S288C, 'SK1'=SK1)
}


normalization_factor <- function(ctrl_input, ctrl_chip,
                                 test_input, test_chip) {
  # Compute Q values
  Q_ctrl_input <- ctrl_input['S288C'] / ctrl_input['SK1']
  Q_ctrl_chip <- ctrl_chip['S288C'] / ctrl_chip['SK1']
  
  Q_test_input <- test_input['S288C'] / test_input['SK1']
  Q_test_chip <- test_chip['S288C'] / test_chip['SK1']
  
  # Compute normalization factors
  a_ctrl <- Q_ctrl_input / Q_ctrl_chip
  a_test <- Q_test_input / Q_test_chip
  
  # Return reference strain-centric normalization factor
  a_test/ a_ctrl
}

#####
setwd('/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1_S288c_Yue_hybrid_MACS2_FE/')
read_counts <- data.frame(
  Condition=c('3h #1', '0h #1'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='Top2-WT0h-747-846-reps_S288C_SK1_Yue_PM_SPMR/stats_HFY77AFXY_n01_AH7797-0h-chipTop2_S288c_SK1_Yue-PM.txt',
         ref_input_counts='Top2-WT0h-747-846-reps_S288C_SK1_Yue_PM_SPMR/stats_HFY77AFXY_n01_AH7797-0h-input_S288c_SK1_Yue-PM.txt',
         test_chip_counts='Top2-WT3h-748-847-reps_S288C_SK1_Yue_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
         test_input_counts='Top2-WT3h-748-847-reps_S288C_SK1_Yue_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt')
  )
)
read_counts
#Condition       NF
#           0h #1 1.000000
#S288C     3h #1 1.641768

#second set
setwd('/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1_S288c_Yue_hybrid_MACS2_FE/')
read_counts <- data.frame(
  Condition=c('3h #2', '0h #2'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='Top2-WT0h-747-846-reps_S288C_SK1_Yue_PM_SPMR/stats_HWV57AFXX_n01_ah7797c0htop2-spikein-0318_S288c_SK1_Yue-PM.txt',
         ref_input_counts='Top2-WT0h-747-846-reps_S288C_SK1_Yue_PM_SPMR/stats_HWV57AFXX_n01_ah7797i0h-0318_S288c_SK1_Yue-PM.txt',
         test_chip_counts='Top2-WT3h-748-847-reps_S288C_SK1_Yue_PM_SPMR/stats_HWV57AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt',
         test_input_counts='Top2-WT3h-748-847-reps_S288C_SK1_Yue_PM_SPMR/stats_HWV57AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt')
  )
)
read_counts
#Condition       NF
#           0h #2 1.000000
#S288C     3h #2 1.279824

nf=(1.641768+1.279824)/2

nfactors=matrix(c(1,1,1.641768,1.279824),4,1)
nfactors=data.frame(nfactors)
nfactors[,2]=c(rep('0h',2),rep('3h',2))
t3_col <- 'black'
t0_col <- 'black'

ggplot(nfactors,aes(x=V2,y=nfactors,fill=V2,width=0.8)) +
  stat_summary(fun.y=mean, geom="bar", width=0.5, alpha=0.25, colour=NA) +
  geom_point(size=1.5, alpha=1) +
  scale_colour_manual('', values=c(t0_col, t3_col,spo11_col), guide=FALSE) +
  scale_fill_manual('', values=c(t0_col, t3_col,spo11_col), guide=FALSE) +
  scale_x_discrete(labels=c(expression('0h'),
                            expression('3h'))) + theme_classic()+
  labs(title = '', x = '', y = 'Relative Top2 amount')
####################################################################################
####################################################################################
# Figure 3b-c
ah7797_3h <- import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1_S288c_Yue_hybrid_MACS2_FE/Top2-WT3h-748-847-reps_S288C_SK1_Yue_PM_SPMR/Top2-WT3h-748-847-reps_S288C_SK1_Yue_PM_SPMR_FE.bdg")
ah7797_0h <- import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1_S288c_Yue_hybrid_MACS2_FE/Top2-WT0h-747-846-reps_S288C_SK1_Yue_PM_SPMR/Top2-WT0h-747-846-reps_S288C_SK1_Yue_PM_SPMR_FE.bdg")

nf_ah7797_3h=(1.641768+1.279824)/2 # calculation in 3a
nf_ah7797_0h=1
normFdiv = function(bdg,nf) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  print(mean(bdg_new$score))
  bdg_new$score <- (bdg_new$score)*nf
  print(mean(bdg_new$score))
  return(bdg_new)
}

ah7797_0h_sk1 <- ah7797_0h[grep("_SK1",ah7797_0h)]
ah7797_3h_sk1 <- ah7797_3h[grep("_SK1",ah7797_3h)]

levelstodrop <- c("chrVI_S288C","chrIII_S288C","chrIV_S288C","chrVIII_S288C","chrII_S288C","chrVII_S288C","chrIX_S288C","chrXVI_S288C","chrXIII_S288C",
                  "chrX_S288C","chrI_S288C","chrV_S288C","chrXV_S288C","chrXII_S288C","chrXIV_S288C","chrXI_S288C")
ah7797_0h_sk1_lev <- dropSeqlevels(ah7797_0h_sk1,levelstodrop)
ah7797_3h_sk1_lev <- dropSeqlevels(ah7797_3h_sk1,levelstodrop)

seqlevels(ah7797_0h_sk1_lev) <- as.character(unlist(strsplit(seqlevels(ah7797_0h_sk1_lev),"_SK1")))
seqlevels(ah7797_3h_sk1_lev) <- as.character(unlist(strsplit(seqlevels(ah7797_3h_sk1_lev),"_SK1")))

ah7797_0hnf = normFdiv(ah7797_0h_sk1_lev,nf_ah7797_0h)
ah7797_3hnf = normFdiv(ah7797_3h_sk1_lev,nf_ah7797_3h)

intergen <- get_intergenic_regions('SK1Yue',as_gr=T)
promoter <- intergen[intergen$type=="divergent"|intergen$type=="tandem"]
mcols(promoter)['widths'] <- width(promoter)
promoter <- promoter[order(promoter$widths,decreasing=T)]
mcols(promoter)['class'] <- DataFrame(class=c(rep(1:3, each=length(promoter)/3),3,3))
midpoint <- floor(width(promoter) / 2)
start(promoter) <- start(promoter) + midpoint
end(promoter) <- start(promoter)

# Figure 3b: 0h
prom1 <- normalizeToMatrix(ah7797_0hnf, promoter[promoter$class==1], value_column = "score",
                           extend = 1000, mean_mode = "weighted", w = 1,empty_value=NA)
prom2 <- normalizeToMatrix(ah7797_0hnf, promoter[promoter$class==2], value_column = "score",
                           extend = 1000, mean_mode = "weighted", w = 1,empty_value=NA)
prom3 <- normalizeToMatrix(ah7797_0hnf, promoter[promoter$class==3], value_column = "score",
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
                     labels = c("-1kb","midpoint","1kb"))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) + geom_line()+ylim(0.8,3.1)
p

# Figure 3c: 3h 
prom13 <- normalizeToMatrix(ah7797_3hnf, promoter[promoter$class==1], value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 1,empty_value=NA)
prom23 <- normalizeToMatrix(ah7797_3hnf, promoter[promoter$class==2], value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 1,empty_value=NA)
prom33 <- normalizeToMatrix(ah7797_3hnf, promoter[promoter$class==3], value_column = "score",
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
                     labels = c("-1kb","midpoint","1kb"))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) + geom_line()+ylim(0.8,3.1)
p

####################################################################################
####################################################################################
# Figure 3d-f

Top2_3 = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1_S288c_Yue_hybrid_MACS2_FE/Top2-WT3h-748-847-reps_S288C_SK1_Yue_PM_SPMR/Top2-WT3h-748-847-reps_S288C_SK1_Yue_PM_SPMR_FE.bdg")
Top2_0 = import_bedGraph("/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveReps_SK1_S288c_Yue_hybrid_MACS2_FE/Top2-WT0h-747-846-reps_S288C_SK1_Yue_PM_SPMR/Top2-WT0h-747-846-reps_S288C_SK1_Yue_PM_SPMR_FE.bdg")
Top2_0_sk1 <- Top2_0[grep("_SK1",Top2_0)]
Top2_3_sk1 <- Top2_3[grep("_SK1",Top2_3)]
levelstodrop <- c("chrVI_S288C","chrIII_S288C","chrIV_S288C","chrVIII_S288C","chrII_S288C","chrVII_S288C","chrIX_S288C","chrXVI_S288C","chrXIII_S288C",
                  "chrX_S288C","chrI_S288C","chrV_S288C","chrXV_S288C","chrXII_S288C","chrXIV_S288C","chrXI_S288C")
Top2_0_sk1_lev <- dropSeqlevels(Top2_0_sk1,levelstodrop)
Top2_3_sk1_lev <- dropSeqlevels(Top2_3_sk1,levelstodrop)
seqlevels(Top2_0_sk1_lev) <- as.character(unlist(strsplit(seqlevels(Top2_0_sk1_lev),"_SK1")))
seqlevels(Top2_3_sk1_lev) <- as.character(unlist(strsplit(seqlevels(Top2_3_sk1_lev),"_SK1")))

normFdiv = function(bdg,nf) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  print(mean(bdg_new$score))
  bdg_new$score <- (bdg_new$score)*nf
  print(mean(bdg_new$score))
  return(bdg_new)
}
nf_ah7797_3h=(1.641768+1.279824)/2 # calculation in 3a
nf_ah7797_0h=1
Top2_3nf = normFdiv(Top2_3_sk1_lev,nf_ah7797_3h)
Top2_0nf = normFdiv(Top2_0_sk1_lev,nf_ah7797_0h)

## Figure 3d: convergent IGRs

intergen <- get_intergenic_regions("SK1Yue",as_gr = T)
conv <- intergen[which(intergen$type=='convergent')]
midpoint <- floor(GenomicRanges::width(conv) / 2)
GenomicRanges::start(conv) <- GenomicRanges::start(conv) + midpoint
GenomicRanges::end(conv) <- GenomicRanges::start(conv)

# find topo signals at hotspots
Top2_3nf_conv <- EnrichedHeatmap::normalizeToMatrix(Top2_3nf, conv,
                                                    extend=1000, w=1,
                                                    mean_mode="weighted",
                                                    value_column="score")
Top2_0nf_conv <- EnrichedHeatmap::normalizeToMatrix(Top2_0nf, conv,
                                                    extend=1000, w=1,
                                                    mean_mode="weighted",
                                                    value_column="score")
Top2_3nf_conv_avrg <- hwglabr2::signal_mean_and_ci(Top2_3nf_conv,
                                                   ci=0.95, rep_bootstrap=1000,
                                                   na_rm=TRUE)
Top2_0nf_conv_avrg <- hwglabr2::signal_mean_and_ci(Top2_0nf_conv,
                                                   ci=0.95, rep_bootstrap=1000,
                                                   na_rm=TRUE)

Top2_3nf_conv_avrg <- data.frame(Data='3h',Position=seq(-999, 1000), Top2_3nf_conv_avrg)
Top2_0nf_conv_avrg <- data.frame(Data='0h',Position=seq(-999, 1000), Top2_0nf_conv_avrg)

# Set up the plot
both <- rbind(Top2_0nf_conv_avrg,Top2_3nf_conv_avrg)
p <- ggplot(both, aes(x=Position, y=Mean, group=Data, fill=Data,color=Data)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  geom_vline(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-1000, 0, 1000),
                     labels = c("-1 kb", "convergent", "1 kb"))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA)+ylim(0.7,3.2)
p + geom_line()

## Figure 3e: promoter-containing IGRs

intergen <- get_intergenic_regions("SK1Yue",as_gr = T)
prom <- intergen[which(intergen$type=='tandem'|intergen$type=='divergent')]
midpoint <- floor(GenomicRanges::width(prom) / 2)
GenomicRanges::start(prom) <- GenomicRanges::start(prom) + midpoint
GenomicRanges::end(prom) <- GenomicRanges::start(prom)

# find topo signals at hotspots
Top2_3nf_prom <- EnrichedHeatmap::normalizeToMatrix(Top2_3nf, prom,
                                                    extend=1000, w=1,
                                                    mean_mode="weighted",
                                                    value_column="score")
Top2_0nf_prom <- EnrichedHeatmap::normalizeToMatrix(Top2_0nf, prom,
                                                    extend=1000, w=1,
                                                    mean_mode="weighted",
                                                    value_column="score")
Top2_3nf_prom_avrg <- hwglabr2::signal_mean_and_ci(Top2_3nf_prom,
                                                   ci=0.95, rep_bootstrap=1000,
                                                   na_rm=TRUE)
Top2_0nf_prom_avrg <- hwglabr2::signal_mean_and_ci(Top2_0nf_prom,
                                                   ci=0.95, rep_bootstrap=1000,
                                                   na_rm=TRUE)

Top2_3nf_prom_avrgdf <- data.frame(Data='3h',Position=seq(-999, 1000), Top2_3nf_prom_avrg)
Top2_0nf_prom_avrgdf <- data.frame(Data='0h',Position=seq(-999, 1000), Top2_0nf_prom_avrg)

# Set up the plot
both <- rbind(Top2_0nf_prom_avrgdf,Top2_3nf_prom_avrgdf)
p <- ggplot(both, aes(x=Position, y=Mean, group=Data, fill=Data,color=Data)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  geom_vline(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-1000, 0, 1000),
                     labels = c("-1 kb", "promoters", "1 kb"))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA)+ylim(0.7,3.2)
p + geom_line()

# Figure 3f: hotspots

SK1Yue_Spo11_DSBs <- get_dsb_hotspots('SK1Yue')
midpoint <- floor(GenomicRanges::width(SK1Yue_Spo11_DSBs) / 2)
GenomicRanges::start(SK1Yue_Spo11_DSBs) <- GenomicRanges::start(SK1Yue_Spo11_DSBs) + midpoint
GenomicRanges::end(SK1Yue_Spo11_DSBs) <- GenomicRanges::start(SK1Yue_Spo11_DSBs)

# find topo signals at hotspots
Top2_3nf_hotspots <- EnrichedHeatmap::normalizeToMatrix(Top2_3nf, SK1Yue_Spo11_DSBs,
                                                        extend=1000, w=1,
                                                        mean_mode="weighted",
                                                        value_column="score")
Top2_0nf_hotspots <- EnrichedHeatmap::normalizeToMatrix(Top2_0nf, SK1Yue_Spo11_DSBs,
                                                        extend=1000, w=1,
                                                        mean_mode="weighted",
                                                        value_column="score")
Top2_3nf_hotspots_avrg <- hwglabr2::signal_mean_and_ci(Top2_3nf_hotspots,
                                                       ci=0.95, rep_bootstrap=1000,
                                                       na_rm=TRUE)
Top2_0nf_hotspots_avrg <- hwglabr2::signal_mean_and_ci(Top2_0nf_hotspots,
                                                       ci=0.95, rep_bootstrap=1000,
                                                       na_rm=TRUE)

Top2_3nf_hotspots_avrg <- data.frame(Data='3h',Position=seq(-999, 1000), Top2_3nf_hotspots_avrg)
Top2_0nf_hotspots_avrg <- data.frame(Data='0h',Position=seq(-999, 1000), Top2_0nf_hotspots_avrg)

# Set up the plot
both <- rbind(Top2_0nf_hotspots_avrg,Top2_3nf_hotspots_avrg)
p <- ggplot(both, aes(x=Position, y=Mean, group=Data, fill=Data,color=Data)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  geom_vline(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-1000, 0, 1000),
                     labels = c("-1 kb", "hotspot", "1 kb"))
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3, color=NA)+ylim(0.7,3.2)
p + geom_line()

