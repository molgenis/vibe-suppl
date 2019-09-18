########
# Name:
# VibeSortingAlgorithmsComparisonResultsGenerator.R
#
# Description:
# Generates plots for comparing the different sorting algorithms within vibe.
# 
# Important:
# Be sure to adjust the items in the config section in "BenchmarkResultsGenerics.R" before running the script!
########

##################
### Libraries  ###
##################

library(rstudioapi)



##################
###    Code    ###
##################

###
### Run R script containing generic R code (only needs to be done once if running multiple plot generating scripts).
###

oldWd <- getwd()
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('./BenchmarkResultsGenerics.R')
setwd(oldWd)



###
### Basic settings & data basic preparation.
###

# Adds R script specific subfolder for storing images.
imgExportDir <- paste0(baseImgExportDir, "vibe_sorting_algorithms/")

# Load data.
vibe.gda_max <- readResultFile("results/vibe_2018-07-06_gda_max.tsv")
vibe.dsi <- readResultFile("results/vibe_2018-07-06_dsi.tsv")
vibe.dpi <- readResultFile("results/vibe_2018-07-06_dpi.tsv")

# Sorts benchmark results so that row order is identical.
vibe.gda_max <- sortRows(vibe.gda_max)
vibe.dsi <- sortRows(vibe.dsi)
vibe.dpi <- sortRows(vibe.dpi)



###
### Calculations for vibe algorithms comparison.
###

# Calculate absolute positions. This is done by processing benchmarkData row-by-row
# and therefore the LOVD row order of positionResults is equal to that of benchmarkData.
positionResults.vibe <- data.frame(vibe.gda_max=resultsPositionCalculator(benchmarkData, vibe.gda_max),
                                   vibe.dsi=resultsPositionCalculator(benchmarkData, vibe.dsi),
                                   vibe.dpi=resultsPositionCalculator(benchmarkData, vibe.dpi))
colnames(positionResults.vibe) <- sapply(sapply(colnames(positionResults.vibe),
                                                strsplit, split=".", fixed=TRUE),
                                         '[', 2)

# Calculate absolute tool rankings.
vibeAlgorithmsRankingCounts <- calculateRankingCounts(positionResults.vibe, calculateRankings(positionResults.vibe))



###
### Plotting figures.
###

toolValuesBoxPlot(positionResults.vibe,
                  'benchmarking_gene_position_vibe_absolute',
                  1, 1000, 500,
                  "",
                  "absolute position of validation genes within the output genes",
                  parMar=c(4.5,5.5,0.5,0.5), ylabMgp=c(4,1,0))

toolRankingBarplot(t(vibeAlgorithmsRankingCounts),
                   'benchmarking_vibe_algorithms_ranking',
                   "", "patient cases", c(0.5,-23))



###
### Removes any script-sepcific variables.
###
rm(imgExportDir, vibe.gda_max, vibe.dsi, vibe.dpi, positionResults.vibe, vibeAlgorithmsRankingCounts)
