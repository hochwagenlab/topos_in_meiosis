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

Top2_level <- data.frame(
  sample=c('Reference','Top2_0h','Top2_0h','Top2_0h','Top2_3h', 'Top2_3h','Top2_3h',"Top2_spo11","Top2_spo11","Top2_rec8","Top2_rec8",'Top2_34c','Top2_34c',"Top2_top2","Top2_top2"),
  Top2levels=c(1,
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2019-11-18_HWC7VAFXY/AH7797-0h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-0h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2019-11-18_HWC7VAFXY/AH7797-0h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i0h-spike_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2018-03-27/AH7797spikein-0h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c0htop2-spikein-0318_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2018-03-27/AH7797spikein-0h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i0h-0318_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2019-01-30_HFY77AFXY/AH7797-0h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-0h-chipTop2_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2019-01-30_HFY77AFXY/AH7797-0h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-0h-input_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2019-11-18_HWC7VAFXY/AH5184-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH5184-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2019-11-18_HWC7VAFXY/AH5184-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH5184-i3h-spike_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2018-03-27/AH5184spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah5184c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2018-03-27/AH5184spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah5184i3h-0318_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2019-11-18_HWC7VAFXY/AH5187-3h_rep1_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH5187-rep1-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2019-11-18_HWC7VAFXY/AH5187-3h_rep1_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH5187-rep1-i3h-spike_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2019-11-18_HWC7VAFXY/AH5187-3h_rep2_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH5187-rep2-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2019-11-18_HWC7VAFXY/AH5187-3h_rep2_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH5187-rep2-i3h-spike_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2019-12-19_H2CVWAFX2/AH7797-34C-3h_spikein-191219_YueSK1_S288C_PM_SPMR/stats_H2CVWAFX2_n01_AH7797-3h-34C-Top2-spike_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2019-12-19_H2CVWAFX2/AH7797-34C-3h_spikein-191219_YueSK1_S288C_PM_SPMR/stats_H2CVWAFX2_n01_AH7797-3h-34C-input-spike_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2019-12-19_H2CVWAFX2/AH7797-34C-3h-rep2_spikein-191219_YueSK1_S288C_PM_SPMR/stats_H2CVWAFX2_n01_AH7797-3h-34C-Top2-spikerep_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2019-12-19_H2CVWAFX2/AH7797-34C-3h-rep2_spikein-191219_YueSK1_S288C_PM_SPMR/stats_H2CVWAFX2_n01_AH7797-3h-34C-input-spikerep_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2019-11-18_HWC7VAFXY/AH7606-3h_rep1_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7606-rep1-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2019-11-18_HWC7VAFXY/AH7606-3h_rep1_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7606-rep1-i3h-spike_S288c_SK1_Yue-PM.txt')
               ),
               suppressMessages(
                 spikein_normalization_factor_from_counts(
                   ref_chip_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-3h-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                                        '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-chipTop2_S288c_SK1_Yue-PM.txt',
                                        '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797c3htop2-spikein-0318_S288c_SK1_Yue-PM.txt'),
                   ref_input_counts=list('2019-11-18_HWC7VAFXY/AH7797-3h_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7797-i3h-spike_S288c_SK1_Yue-PM.txt',
                                         '2019-01-30_HFY77AFXY/AH7797-3h_spikein-013019_YueSK1_S288C_PM_SPMR/stats_HFY77AFXY_n01_AH7797-3h-inputTop2_S288c_SK1_Yue-PM.txt',
                                         '2018-03-27/AH7797spikein-3h-032718_YueSK1_S288C_PM_SPMR_PE/stats_HW5G7AFXX_n01_ah7797i3h-0318_S288c_SK1_Yue-PM.txt'),
                   test_chip_counts='2019-11-18_HWC7VAFXY/AH7606-3h_rep2_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7606-rep2-Top2chip-spike_S288c_SK1_Yue-PM.txt',
                   test_input_counts='2019-11-18_HWC7VAFXY/AH7606-3h_rep2_spikein-191118_YueSK1_S288C_PM_SPMR/stats_HWC7VAFXY_n01_AH7606-rep2-i3h-spike_S288c_SK1_Yue-PM.txt')
               )))


ggplot(Top2_level[which(Top2_level$sample!='Reference'),],aes(x=sample,y=Top2levels,fill=sample,width=0.8)) +
  stat_summary(fun.y=mean, geom="bar", width=0.5, alpha=0.25, colour=NA) +
  geom_point(size=1.5, alpha=1) +
  scale_x_discrete(labels=c(expression('0h'),
                            expression('3h'),expression('34C'),expression('spo11'),expression('rec8'),expression('top2-1'))) +
  theme_classic()  + theme(legend.position = "none") +
  labs(title = '', x = '', y = 'Relative Top2 amount')

####################################################################################
####################################################################################
# Figure 3b

spo11oligo <- rtracklayer::import.bedGraph("Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg")
gff <- hwglabr2::get_gff('SK1Yue')
transcription <- read.csv('2016.03.16-2h+3h/2017.06.16_SK1Yue_EdgeR_tpm.csv')

gff <- gff[which(gff$type=='gene')]
colnames(transcription)[1] <- "ID"
gff <- data.frame(gff)
gff_txn <- merge(x=gff,y=transcription[,c(1,3)],by='ID', all.x = TRUE)
gff_txn <- gff_txn[which(gff_txn$seqnames!='chrMT'),]
gff_txn <- gff_txn[which(gff_txn$seqnames!='scplasm1'),]

gffgr <- GRanges(seqnames = gff_txn$seqnames,IRanges(start = gff_txn$start,end=gff_txn$end),strand = gff_txn$strand,transcription=gff_txn$AH119_3h,gene=gff_txn$Name)
negstrand=gffgr[strand(gffgr)=='-']
posstrand=gffgr[strand(gffgr)=='+']
negstrandswitch <- GRanges(seqnames = seqnames(negstrand),ranges = IRanges(start=end(negstrand),end=end(negstrand)),strand=strand(negstrand),transcription=negstrand$transcription,gene=negstrand$gene)
end(posstrand) <- start(posstrand)
gffall <- c(posstrand,negstrandswitch)
gffall_sort <- gffall[order(gffall$transcription,decreasing = T)]
gffall_sort <- gffall_sort[which(gffall_sort$transcription!='<NA>')]
mcols(gffall_sort) <- DataFrame(class=c(rep(1:4, each=length(gffall_sort)/4),4,4))

spo11oligo_atg <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, gffall_sort,
                                                     extend=c(500,500), w=1,empty_value=NA,
                                                     mean_mode="weighted",
                                                     value_column='score')
col_fun <- colorRamp2(quantile(spo11oligo_atg, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
partition <- gffall_sort$class
EnrichedHeatmap(spo11oligo_atg, col = col_fun, name = "Spo11", row_title_rot = 0,
                row_order = 1:length(gffall_sort),
                split=gffall_sort$class,
                axis_name = c("-500", "ATG","500"))+
  Heatmap(partition, col = structure(1:4, names = as.character(1:4)), name = "",row_order = 1:length(gffall_sort),
          show_row_names = FALSE, width = unit(5, "mm"))

# average lines
spo11_1 <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, gffall_sort[which(gffall_sort$class == 1)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
spo11_1ci <- hwglabr2::signal_mean_and_ci(signal_data=spo11_1,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
spo11_2 <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, gffall_sort[which(gffall_sort$class == 2)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
spo11_2ci <- hwglabr2::signal_mean_and_ci(signal_data=spo11_2,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
spo11_3 <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, gffall_sort[which(gffall_sort$class == 3)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
spo11_3ci <- hwglabr2::signal_mean_and_ci(signal_data=spo11_3,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
spo11_4 <- EnrichedHeatmap::normalizeToMatrix(spo11oligo, gffall_sort[which(gffall_sort$class == 4)],
                                              extend=1000, w=1,
                                              mean_mode="weighted",
                                              value_column="score")
spo11_4ci <- hwglabr2::signal_mean_and_ci(signal_data=spo11_4,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
group1_gg <- data.frame(Data="1",Position=seq(1, 2000), spo11_1ci)
group2_gg <- data.frame(Data="2",Position=seq(1, 2000), spo11_2ci)
group3_gg <- data.frame(Data="3",Position=seq(1, 2000), spo11_3ci)
group4_gg <- data.frame(Data="4",Position=seq(1, 2000), spo11_4ci)
allgroups <- rbind(group1_gg,group2_gg,group3_gg,group4_gg)

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
# Figure 3c

# Spo11 oligos
spo11oligo <- rtracklayer::import.bedGraph("Spo11oligo_WT1_SRR-clip-MACS2_extsize37_treat_pileup.bdg")

intergen <- hwglabr2::get_intergenic_regions("SK1Yue",as_gr=T)
prom <- intergen[intergen$type=="divergent"|intergen$type=='tandem']
prom <- prom[order(width(prom),decreasing = T)]
midpoint <- floor(width(prom) / 2)
start(prom) <- start(prom) + midpoint
end(prom) <- start(prom)

prommat <- normalizeToMatrix(spo11oligo, prom, value_column = "score",
                             extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
col_fun <- colorRamp2(quantile(prommat, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
EnrichedHeatmap(prommat, col = col_fun, name = "Spo11",
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = "purple4"),
                                                                         show_error = TRUE)),
                top_annotation_height = unit(2, "cm"), row_title_rot = 0,
                axis_name = c("-1 kb", "promoters", "1 kb"),
                row_order = 1:length(prom))

# Nucleosomes
mnase3 = rtracklayer::import.bedGraph("Nucleosome_reps-SK1-MACS2_treat_pileup.bdg")
gendiv = function(bdg) {
  gavg = average_chr_signal(bdg)$genome_avrg
  print(gavg)
  bdg_new <- bdg
  bdg_new$score <- bdg_new$score/gavg
  return(bdg_new)
}
mnase3d = gendiv(mnase3)
intergen <- hwglabr2::get_intergenic_regions("SK1Yue",as_gr=T)
prom <- intergen[intergen$type=="divergent"|intergen$type=='tandem']
prom <- prom[order(width(prom),decreasing = T)]
midpoint <- floor(width(prom) / 2)
start(prom) <- start(prom) + midpoint
end(prom) <- start(prom)

prommat <- normalizeToMatrix(mnase3d, prom, value_column = "score",
                             extend = 1000, mean_mode = "weighted", w = 10,empty_value=NA)
col_fun <- colorRamp2(quantile(prommat, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
EnrichedHeatmap(prommat, col = col_fun, name = "Nuc 3h",
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = "purple4"),
                                                                         show_error = TRUE)),
                top_annotation_height = unit(2, "cm"), row_title_rot = 0,
                axis_name = c("-1 kb", "promoters", "1 kb"),
                row_order = 1:length(prom))

####################################################################################
####################################################################################
# Figure 3d

extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "numeric")
import.broadPeak <- function(...) {
  rtracklayer::import(..., format="BED", extraCols=extraCols_broadPeak)
}

SK1Yue_Spo11_DSB <- get_dsb_hotspots('SK1Yue')
Top2_peak = import.broadPeak("Top2-wildtype-413-504-Reps-SK1Yue-PM_B3W3_MACS2_peaks.broadPeak")
Top1_peak = import.broadPeak("AH9847Myc-3h-735-841-Reps-SK1Yue-PM_B3W4_MACS2_peaks.broadPeak")

#Hotspots that overlap with Top2
Spo11Top2 = findOverlaps(Top2_peak,SK1Yue_Spo11_DSB)
Spo11Top2 = SK1Yue_Spo11_DSB[subjectHits(Spo11Top2)]
#Hotspots that overlap with Top1
Spo11Top1 = findOverlaps(Top1_peak,SK1Yue_Spo11_DSB)
Spo11Top1 = SK1Yue_Spo11_DSB[subjectHits(Spo11Top1)]

#Find Top2 hotspots that overlap with Top1
Spo11Top2Top1 = findOverlaps(Spo11Top1,Spo11Top2)
Spo11Top2Top1 = Spo11Top2[subjectHits(Spo11Top2Top1)] #n=1694

#find Top2 peaks that do not overlap with the rest
Top2alone = Spo11Top2[!(Spo11Top2 %over% Spo11Top1)] #n=780

#find Top1 peaks that do not overlap with the rest
Top1alone = Spo11Top1[!(Spo11Top1 %over% Spo11Top2)] #n=253

#hotspots that do not overlap with Top1 or Top2
SK1Yue_Spo11_DSBalone = SK1Yue_Spo11_DSB[!(SK1Yue_Spo11_DSB %over% Spo11Top1)]
SK1Yue_Spo11_DSBalone = SK1Yue_Spo11_DSBalone[!(SK1Yue_Spo11_DSBalone %over% Spo11Top2)] # n=464

alldata=list(SK1Yue_Spo11_DSBalone$score, #no overlaps
             Top1alone$score, #Top1 alone
             Top2alone$score, #Top2 alone
             Spo11Top2Top1$score) #Top1+Top2
boxplot(alldata,names=c("None", "Top1 Only", "Top2 Only","Top1&Top2"),outline=F,ylab="Signal Intensity",frame.plot=F,cex.lab=1.5,cex.axis=1.25)

p <- c(wilcox.test(SK1Yue_Spo11_DSBalone$score,Top1alone$score,paired = F)$p.value,
       wilcox.test(SK1Yue_Spo11_DSBalone$score,Top2alone$score,paired = F)$p.value,
       wilcox.test(SK1Yue_Spo11_DSBalone$score,Spo11Top2Top1$score,paired = F)$p.value,
       wilcox.test(Top1alone$score,Top2alone$score,paired = F)$p.value,
       wilcox.test(Top1alone$score,Spo11Top2Top1$score,paired = F)$p.value,
       wilcox.test(Top2alone$score,Spo11Top2Top1$score,paired = F)$p.value)
p.adjust(p, method = "bonferroni", n = length(p))
