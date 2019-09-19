########
# Name:
# VibeAdjustmentsComparison.R
#
# Description:
# Generates plots for comparing the different tools.
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
imgExportDir <- paste0(baseImgExportDir, "vibe_adjustments/")

# Defines color palette.
adjustmentColors <- brewer.pal(4, 'Set2') # Adjust if changing number of tools!!!

# Load data.
vibe.v5.simple <- readResultFile("results/vibe_2019-09-15_none.tsv")
vibe.v5.pda <- readResultFile("results/vibe_2018-07-06_none.tsv")
vibe.v6.simple <- readResultFile("results/vibe_2019-09-19_v6-simple_none.tsv")
vibe.v6.pda <- readResultFile("results/vibe_2019-09-19_v6-old_none.tsv")

# Sorts benchmark results so that row order is identical.
vibe.v5.simple <- sortRows(vibe.v5.simple)
vibe.v5.pda <- sortRows(vibe.v5.pda)
vibe.v6.simple <- sortRows(vibe.v6.simple)
vibe.v6.pda <- sortRows(vibe.v6.pda)



###
### Calculations for comparison of vibe changes - positions
###

# Calculate absolute positions. This is done by processing benchmarkData row-by-row
# and therefore the LOVD row order of positionResults is equal to that of benchmarkData.
positionResults <- data.frame(vibe.v5.simple=resultsPositionCalculator(benchmarkData, vibe.v5.simple),
                              vibe.v5.pda=resultsPositionCalculator(benchmarkData, vibe.v5.pda),
                              vibe.v6.simple=resultsPositionCalculator(benchmarkData, vibe.v6.simple),
                              vibe.v6.pda=resultsPositionCalculator(benchmarkData, vibe.v6.pda))

# Calculate total number of genes.
totalResults <- data.frame(vibe.v5.simple=calculateTotalGenesFound(vibe.v5.simple),
                           vibe.v5.pda=calculateTotalGenesFound(vibe.v5.pda),
                           vibe.v6.simple=calculateTotalGenesFound(vibe.v6.simple),
                           vibe.v6.pda=calculateTotalGenesFound(vibe.v6.pda),
                           row.names=rownames(vibe.v5.simple))

# Replicates some of the totalResults so that size is equal to positionResults.
totalResults <- totalResults[benchmarkData$lovd,]


# Calculate relative positions.
relativePositionResults <- positionResults / totalResults



###
### Calculations for comparison of different tools (with best scoring vibe result) - found matches within position cutoffs
###

# Calculates how many of the genes are found when only looking within specific
# absolute cutoffs of the total number of hits available.
highestPossibleResult <- apply(totalResults, 2, max)
cutoffRanges <- seq(1, max(highestPossibleResult), 1) # how specific looking for cutoffs
genesFoundWithinAbsoluteCutoff <- apply(positionResults, 2, function(toolAbsoluteScores, range=cutoffRanges) {
  sapply(range, function(cutoff, scores=toolAbsoluteScores) {
    length(which(scores < cutoff))
  })
})
rownames(genesFoundWithinAbsoluteCutoff) <- cutoffRanges

# Replaces values with NA for positions if higher than highest available
# totalResult for that tool.
for(x in names(highestPossibleResult)) {
  if(highestPossibleResult[x] < nrow(genesFoundWithinAbsoluteCutoff)) {
    genesFoundWithinAbsoluteCutoff[(highestPossibleResult[x]+1):nrow(genesFoundWithinAbsoluteCutoff),x] <- NA
  }
}
rm(x, highestPossibleResult)

# Calculates how many of the genes are found when only looking within specific
# relative cutoffs of the total number of hits available.
cutoffRanges <- seq(0.0001, 1, 0.0001) # how specific looking for cutoffs
genesFoundWithinRelativeCutoff <- apply(relativePositionResults, 2, function(toolRelativeScores, range=cutoffRanges) {
  sapply(range, function(cutoff, scores=toolRelativeScores) {
    length(which(scores < cutoff))
  })
})
rownames(genesFoundWithinRelativeCutoff) <- cutoffRanges

# Converts above two into percentage of found hits compared to total hits.
genePercentageFoundWithinAbsoluteCutoff <- genesFoundWithinAbsoluteCutoff / nrow(positionResults)
genePercentageFoundWithinRelativeCutoff <- genesFoundWithinRelativeCutoff / nrow(relativePositionResults)

# Removes variable so that it cannot accidentally be used later on.
rm(cutoffRanges)



###
### Curves showing how many genes were found within different cutoffs.
###

# Plots how often genes were found using cutoffs of total available hits.
# Shows how many hits were found when looking at a specific absolute cutoff
# max positions from all genes found.
plotMatchesFoundWihinRangeCutoff(genePercentageFoundWithinAbsoluteCutoff[1:5000,],
                                 'vibe_found_genes_for_absolute_cutoffs',
                                 "absolute cutoff from the output genes",
                                 "fraction of patient cases for which the gene was found",
                                 adjustmentColors, 500, addAbline=FALSE)

# Shows how many hits were found when looking at a specific cutoff fractions
# from all genes found.
plotMatchesFoundWihinRangeCutoff(genePercentageFoundWithinRelativeCutoff,
                                 'vibe_found_genes_for_relative_cutoffs',
                                 "relative cutoff from the output genes",
                                 "fraction of patient cases for which the gene was found",
                                 adjustmentColors, 0.1)



###
### Removes any script-sepcific variables.
###
rm(adjustmentColors, imgExportDir, vibe.v5.simple, vibe.v5.pda,
   genePercentageFoundWithinAbsoluteCutoff, genePercentageFoundWithinRelativeCutoff,
   genesFoundWithinAbsoluteCutoff, genesFoundWithinRelativeCutoff,
   positionResults, relativePositionResults, totalResults)
