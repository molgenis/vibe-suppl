########
# Name:
# VibeVersionsComparison
#
# Description:
# Generates plots for comparing the different versions of vibe.
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
imgExportDir <- paste0(baseImgExportDir, "vibe_versions/")

# Defines color palette.
#versionColors <- brewer.pal(5, 'Set2') # Adjust if changing number of tools!!!
versionColors <- c('slategray2','steelblue4')

# Load data and sorts it (ensures identical row order).
vibe.1_0.results <- sortRows(readResultFile("vibe_50ecddcff7436a5023fa4fbe65098c45fa76d2ce_results.tsv"))
vibe.1_0.time <- sortRows(readTimesFile("vibe_50ecddcff7436a5023fa4fbe65098c45fa76d2ce_time.tsv"))

vibe.2_0_results <- sortRows(readResultFile("vibe_934b26a5c8d12fbe36e8ef63da945eae21217bfb_results_converted.tsv"))
vibe.2_0.time <- sortRows(readTimesFile("vibe_934b26a5c8d12fbe36e8ef63da945eae21217bfb_time.tsv"))

###
### Adjustments/calculations for times.
###

# Merges times into dataframe.
vibe_times <- data.frame(old=vibe.1_0.time[[1]],
                         new=vibe.2_0.time[[1]],
                         row.names=rownames(vibe.1_0.time))

timeGainFactors <- vibe_times$old / vibe_times$new

print(round(min(timeGainFactors), 2))
print(round(max(timeGainFactors), 2))
print(round(mean(timeGainFactors), 2))


###
### Calculations for comparison of vibe changes - positions
###

# Calculate absolute positions. This is done by processing benchmarkData row-by-row
# and therefore the LOVD row order of positionResults is equal to that of benchmarkData.
positionResults <- data.frame(old=resultsPositionCalculator(benchmarkData, vibe.1_0.results),
                              new=resultsPositionCalculator(benchmarkData, vibe.2_0_results))

# Calculate total number of genes.
totalResults <- data.frame(old=calculateTotalGenesFound(vibe.1_0.results),
                           new=calculateTotalGenesFound(vibe.2_0_results),
                           row.names=rownames(vibe.1_0.results))

# Replicates some of the totalResults so that size is equal to positionResults.
totalResults <- totalResults[benchmarkData$lovd,]

# Calculate relative positions.
relativePositionResults <- positionResults / totalResults



###
### Calculations for comparison of tool adjustments - found matches within position cutoffs
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
### Shows cases where at least 1 of the versions has NA.
###

positionResults[!complete.cases(positionResults),]



###
### Plots differences in time.
###
dataToPlot <- vibe_times[order(vibe_times$old, decreasing=TRUE),]

initializeGraphicsDevice("vibe_time", width=10, height=5)
barplot(dataToPlot$old, ylim=c(0,25),
        border=NA, col=versionColors[1], space=0, las=1,
        xlab='the individual benchmarks (sorted on time needed by old version)',
        ylab='user+sys time (minutes)')
barplot(dataToPlot$new, ylim=c(0,25), border=NA,
        col=versionColors[2], space=0, las=1, add=T, yaxt='n')
legend(nrow(dataToPlot)*0.8,25, c("old","new"), fill=versionColors, bty = "n")
dev.off()

rm(dataToPlot)



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
                                 versionColors, 500, addAbline=FALSE)

# Shows how many hits were found when looking at a specific cutoff fractions
# from all genes found.
plotMatchesFoundWihinRangeCutoff(genePercentageFoundWithinRelativeCutoff,
                                 'vibe_found_genes_for_relative_cutoffs',
                                 "relative cutoff from the output genes",
                                 "fraction of patient cases for which the gene was found",
                                 versionColors, 0.1)



###
### Compares individual items between old and new version.
###

# Plots vibe against amelie (absolute).
dataToPlot <- log10(data.frame(old=positionResults$old,
                               new=positionResults$new))
plotToolComaprison(dataToPlot,
                   'benchmark_positions_old_new_absolute',
                   max(dataToPlot, na.rm=TRUE),
                   "absolute position of validation genes old version (log10)",
                   "absolute position of validation genes new version (log10)")

# Plots vibe against amelie (relative).
dataToPlot <- log10(data.frame(old=relativePositionResults$old,
                               new=relativePositionResults$new))
plotToolComaprison(dataToPlot,
                   'benchmark_positions_old_new_relative',
                   xyMin=floor(min(dataToPlot, na.rm=TRUE)),
                   xyMax=ceiling(max(dataToPlot, na.rm=TRUE)),
                   "relative position of validation genes old version (log10)",
                   "relative position of validation genes new version (log10)")
