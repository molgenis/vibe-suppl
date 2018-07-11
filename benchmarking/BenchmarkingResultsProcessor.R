########
# Name:
# BenchmarkingResultsProcessor.R
#
# Description:
# 
########

##################
### Libraries  ###
##################

#install.packages("VennDiagram")
library(VennDiagram)
library(RColorBrewer)





##################
### Functions  ###
##################

########
# Name:
# readResultFile
#
# Description:
# Reads in a single file containing benchmarking results.
#
# Input:
# filePath -  The path to the file to be loaded.
#
# Output:
# A table from the result data. Though as the input is expected to only contain
# 2 columns and one of these is used as row names, type becomes a list.
########
readResultFile <- function(filePath) {
  read.table(paste0(baseDir, filePath),
             header=T, sep="\t", colClasses=c("character"),
             row.names=1)
}

########
# Name:
# resultsPositionCalculator
#
# Description:
# For each row in benchmarkData (containing an LOVD and a gene), looks at the
# benchmarkResults to determine the gene positions.
#
# Input:
# benchmarkData - The data on which the benchmark is based upon. Should at least
#                 contain a column "gene" containing the gene to be found
#                 (determening the gene position) and "lovd" to determine the ID
#                 of which phenotype-set was used.
# benchmarkResults - All the output from a benchmark. Each row contains the LOVD
#                    as row name and the ordered genes as a column named "genes"
#                    (these genes are comma-separated within this single field).
#
# Output:
# 
########
resultsPositionCalculator <- function(benchmarkData, benchmarkResults) {
  apply(benchmarkData, 1, singleResultPositionCalculator,
        benchmarkResults=benchmarkResults)
}

########
# Name:
# singleResultPositionCalculator
#
# Description:
# Processes a single row of the benchmarkData as defined in the function
# resultsPositionCalculator. For this single benchmark slice, it looks up the
# corresponding LOVD in benchmarkResults and then defines the gene position.
#
# Input:
# benchmarkDataRow - A list (row from benchmarkData given by
#                    resultsPositionCalculator) containing an LOVD and a gene.
# benchmarkResults - All the output from a benchmark. Each row contains the LOVD
#                    as row name and the ordered genes as a column named "genes"
#                    (these genes are comma-separated within this single field).
#
# Output:
# 
########
singleResultPositionCalculator <- function(benchmarkDataRow, benchmarkResults) {
  match(benchmarkDataRow["gene"],
        strsplit(benchmarkResults[benchmarkDataRow["lovd"],"suggested_genes"], ",")[[1]])
}

########
# Name:
# calculateTotalGenesFound
#
# Description:
# Calculates the total number of genes per item for the input list benchmarkResults.
#
# Input:
# benchmarkResults - All the output from a benchmark. Each row contains the LOVD
#                    as row name and the ordered genes as a column named "genes"
#                    (these genes are comma-separated within this single field).
#
# Output:
# An integer vector containing the total number of genes for each (named) item
# (in the row order of benchmarkResults).
########
calculateTotalGenesFound <- function(benchmarkResults) {
  unlist(lapply(sapply(benchmarkResults[,"suggested_genes"], strsplit, split=","), length), use.names=FALSE)
}

########
# Name:
# sortRows
#
# Description:
# Orders the rows based on their row name.
#
# Input:
# benchmarkResults - All the output from a benchmark. Each row contains the LOVD
#                    as row name and the ordered genes as a column named "genes"
#                    (these genes are comma-separated within this single field).
#
# Output:
# benchmarkResults with the rows ordered on their name.
########
sortRows <- function(benchmarkResults) {
  benchmarkResults[order(as.numeric(rownames(benchmarkResults))), , drop=FALSE]
}

########
# Name:
# calculateRankings
#
# Description:
# 
#
# Input:
# 
#
# Output:
#
########
calculateRankings <- function(data) {
  positionNames <- c("first", "second", "third", "fourth", "fifth")
  
  # Generates character matrix with on the rows the different benchmark input sets
  # and on the columns the position of the tool (1st column is best ranking tool,
  # second column is second ranking tool, et cetera).
  rankingResults <- sapply(apply(data, 1, sort, na.last=NA), names)
  rankingResults <- sapply(1:ncol(data), function(x) { sapply(rankingResults, '[', x) })
  colnames(rankingResults) <- positionNames[1:ncol(data)]
  
  return(rankingResults)
}

########
# Name:
# calculateRankingCounts
#
# Description:
# 
#
# Input:
# 
#
# Output:
#
########
calculateRankingCounts <- function(data, rankingResults) {
  # Counts how often each tools scores a specific position (first, second, et cetera).
  rankingResultCounts <- apply(rankingResults, 2, function(x) {
    table(factor(x, levels=colnames(data))) # Factor ensures zeros can be "counted" with table.
  })
  
  # Count missing values.
  resultsWithNa <- data[apply(apply(data, c(1,2), is.na), 1, any),]
  naCounts <- apply(apply(resultsWithNa, c(1,2), is.na), 2, sum)
  
  # Generates complete table for further usage.
  rankingResultCounts <- cbind(rankingResultCounts, "no hit"=naCounts)
  
  return(rankingResultCounts)
}

########
# Name:
# toolValuesBoxPlot
#
# Description:
# 
#
# Input:
# 
#
# Output:
#
########
toolValuesBoxPlot <- function(data, file, firstYValue, yAxisFreq, yAxisAblineFreq, main, ylab) {
  postscript(file, width=7, height=8)
  
  yAxisMax <- ceiling(max(data, na.rm=T)/yAxisFreq)*yAxisFreq
  boxplot(data, xaxt='n',yaxt='n', pch="")
  abline(h=firstYValue, col="gray92")
  abline(h=seq(yAxisAblineFreq, yAxisMax, yAxisAblineFreq), col="gray92")
  boxplot(data, las=1, pch=20, yaxt='n',
          main=main, xlab="tool", ylab=ylab, col="white", add=T)
  axis(2, las=1, at=firstYValue)
  axis(2, las=1, at=seq(yAxisFreq, yAxisMax,yAxisFreq))
  
  dev.off()
}

########
# Name:
# toolRankingBarplot
#
# Description:
# 
#
# Input:
# 
#
# Output:
#
########
toolRankingBarplot <- function(data, file, main, legendPos) {
  GreenToRedColors <- rev(brewer.pal(nrow(data), 'RdYlGn'))
  postscript(file, width=7, height=4)
  
  barplot(data, col=GreenToRedColors, las=1,
          main=main)
  par(xpd=TRUE) # no clipping for drawing outside plot
  legend(legendPos[1], legendPos[2], rownames(data),
         fill=GreenToRedColors, ncol=nrow(data))
  
  dev.off()
}

########
# Name:
# benchmarkQuintupleVenn
#
# Description:
# Draws a venn diagram using data as input and writes text outside of it
# indicating the number from outside to show how many were not included by any
# of the items from the venn diagram.
#
# Input:
# data - A list of vectors (e.g., integers, chars), with each component
#        corresponding to a separate circle in the Venn diagram (direct citation
#        from venn.diagram() from the library VennDiagram).
# outside - A number indicating how many hits were outside the venn diagram.
# file - Full path (including the ".eps" extension) to the location where the
#        plot should be stored on local storage.
# colors - The colors to be used for the venn diagram.
#
# Output:
#
########
benchmarkQuintupleVenn <- function(data, outside, file, colors) {
  cairo_ps(file, width=7, height=7)
  
  grid.draw(venn.diagram(data, NULL,
                         fontfamily="Helvetica", main.fontfamily="Helvetica",
                         sub.fontfamily="Helvetica", cat.fontfamily="Helvetica",
                         fill=colors, cat.just=list(c(0.5,1),
                                                    c(-0.5,-6),
                                                    c(0,0),
                                                    c(1,1),
                                                    c(1.5,-6))
  ))
  grid.text(paste("none\n", outside), 0.1, 0.2)
  
  dev.off()
}





##################
###    Code    ###
##################

###
### Basic settings & data basic preparation.
###

# Defaults
oldPar <- par()
setEPS() # Sets EPS engine for writing images.

# Adjust to local system.
baseDir <- '~/programming/projects/vibe/data/benchmarking/'
imgExportDir <- '~/Documents/afstuderen/verslag/img/'

# Load data.
benchmarkData <- read.table(paste0(baseDir,"benchmark_data.tsv"), header=T,
                            sep="\t",colClasses=c(rep("character", 3),
                                                  "factor", "character"))

amelie <- readResultFile("results/amelie.tsv")
geneNetwork <- readResultFile("results/gene_network.tsv")
phenomizer <- readResultFile("results/phenomizer.tsv")
phenotips <- readResultFile("results/phenotips.tsv")
vibe.gda_max <- readResultFile("results/vibe_2018-07-06_gda_max.tsv")
vibe.dsi <- readResultFile("results/vibe_2018-07-06_dsi.tsv")
vibe.dpi <- readResultFile("results/vibe_2018-07-06_dpi.tsv")

# Sorts benchmark results so that row order is identical.
amelie <- sortRows(amelie)
geneNetwork <- sortRows(geneNetwork)
phenomizer <- sortRows(phenomizer)
phenotips <- sortRows(phenotips)
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
### Plotting vibe-only figures.
###

toolValuesBoxPlot(positionResults.vibe,
                  paste0(imgExportDir, 'benchmarking_gene_position_vibe_absolute.eps'),
                  1, 500, 500,
                  "position of relevant genes among different vibe sorting algorithms",
                  "position")

toolRankingBarplot(t(vibeAlgorithmsRankingCounts),
                   paste0(imgExportDir, 'benchmarking_vibe_algorithms_ranking.eps'),
                   "gene position ranked among the tools", c(0.5,-80))



###
### Calculations for comparison of different tools (with best scoring vibe result).
###

# Sets which vibe algorithm should be used for further processing.
vibe <- vibe.gda_max


# Calculate absolute positions. This is done by processing benchmarkData row-by-row
# and therefore the LOVD row order of positionResults is equal to that of benchmarkData.
positionResults <- data.frame(amelie=resultsPositionCalculator(benchmarkData, amelie),
                              geneNetwork=resultsPositionCalculator(benchmarkData, geneNetwork),
                              phenomizer=resultsPositionCalculator(benchmarkData, phenomizer),
                              phenotips=resultsPositionCalculator(benchmarkData, phenotips),
                              vibe=resultsPositionCalculator(benchmarkData, vibe))

# Calculate total number of genes.
totalResults <- data.frame(amelie=calculateTotalGenesFound(amelie),
                           geneNetwork=calculateTotalGenesFound(geneNetwork),
                           phenomizer=calculateTotalGenesFound(phenomizer),
                           phenotips=calculateTotalGenesFound(phenotips),
                           vibe=calculateTotalGenesFound(vibe),
                           row.names=rownames(amelie))

# Replicates some of the totalResults so that size is euql to positionResults.
totalResults <- totalResults[benchmarkData$lovd,]


# Calculate relative positions.
relativePositionResults <- positionResults / totalResults

# Calculate absolute tool rankings.
toolRankings <- calculateRankings(positionResults)
toolRankingCounts <- calculateRankingCounts(positionResults, toolRankings)

# Retrieves all used input phenotypes.
allPhenotypeNames <- sort(unique(unlist(sapply(benchmarkData[,5], strsplit, ";"), use.name=FALSE)))

# Generates a dataframe with 0-values to be filled in later for the phenotypes.
# This specific dataframe is focused on counting the phenotypes for when a tool
# ranked first among all tools.
phenotypeCountsWhenToolRankedFirst <- data.frame(amelie=rep(0,length(allPhenotypeNames)),
                                                 geneNetwork=rep(0,length(allPhenotypeNames)),
                                                 phenomizer=rep(0,length(allPhenotypeNames)),
                                                 phenotips=rep(0,length(allPhenotypeNames)),
                                                 vibe.20180503=rep(0,length(allPhenotypeNames)),
                                                 row.names=allPhenotypeNames)

# Looks per tool for the phenotype-sets in which it ranked the best and from these
# merges all phenotypes and looks at how often each phenotypes occur.
foundCounts <-
  sapply(rownames(toolRankingCounts), function(toolName) {
    # All phenotype frequencies for a single tool in the cases it ranked best.
    table(unlist(
      # toolRankingResults[,"first"]==toolName filters benchmarkData on the
      # phenotype-sets where toolName ranked "first".
      # 5 is the phenotype names column in benchmarkData.
      sapply(benchmarkData[which(toolRankings[,"first"]==toolName),5],
             strsplit, ";"),
      use.names=FALSE))
  })

# Adds the results from foundCounts to phenotypeCountsWhenToolRankedFirst.
for(toolName in colnames(phenotypeCountsWhenToolRankedFirst)) {
  phenotypeCountsWhenToolRankedFirst[names(foundCounts[[toolName]]),toolName] <- foundCounts[[toolName]]
}
rm(foundCounts)

# The frequencies of the input phenotypes for which a tool ranked first (no all NA).
phenotypeTotals <- apply(phenotypeCountsWhenToolRankedFirst, 1, sum)



###
### Data info.
###

# Check whether the rows of all benchmarks are identically ordered.
rownames(amelie) == rownames(geneNetwork) &&
  rownames(geneNetwork) == rownames(phenomizer) &&
  rownames(phenomizer) == rownames(phenotips) &&
  rownames(phenotips) == rownames(vibe)

# Check whether there is any phenotype set for which no single gene was found.
any(is.na(totalResults))

# Check whether tool ranking return values indicating they are complete (and equal!).
apply(toolRankingCounts, 1, sum)



###
### Plotting figures.
###

# Generics.
toolColors <- brewer.pal(ncol(positionResults), 'Set3')

# Boxplot comparing absolute positions of genes.
toolValuesBoxPlot(positionResults,
                  paste0(imgExportDir, 'benchmarking_gene_position_absolute.eps'),
                  1, 500, 100,
                  "position of relevant genes among different tools",
                  "position")

# Boxplot comparing number of suggested genes found.
toolValuesBoxPlot(totalResults,
                  paste0(imgExportDir, 'benchmarking_gene_position_total_hits.eps'),
                  1, 2000, 500,
                  "the number of output genes among different tools",
                  "number of genes")

# Boxplot comparing relative positions of genes (absolute position / number of suggested genes found).
toolValuesBoxPlot(relativePositionResults,
                  paste0(imgExportDir, 'benchmarking_gene_position_relative.eps'),
                  0.0, 0.2, 0.05,
                  "relative position of relevant genes among different tools",
                  "relative position")

# Barplot showing the tool rankings.
toolRankingBarplot(t(toolRankingCounts),
                   paste0(imgExportDir, 'benchmarking_tool_ranking.eps'),
                   "gene position ranked among the tools", c(0,-80))

# Plot differences in whether the gene was found.
benchmarkQuintupleVenn(apply(!is.na(positionResults), 2, which),
                       sum(apply(is.na(positionResults), 1, all)),
                       paste0(imgExportDir, 'benchmarking_overlap_genes_found_absolute.eps'),
                       toolColors)

# Plot differences in whether the gene was found within the first 100 positions.
benchmarkQuintupleVenn(apply(positionResults <=1000, 2, which),
                       sum(apply(positionResults > 1000, 1, all, na.rm=T)),
                       paste0(imgExportDir, 'benchmarking_overlap_genes_found_absolute_max_1000.eps'),
                       toolColors)

# Plot differences in whether the gene was found within the first 100 positions.
benchmarkQuintupleVenn(apply(positionResults <=100, 2, which),
                       sum(apply(positionResults > 100, 1, all, na.rm=T)),
                       paste0(imgExportDir, 'benchmarking_overlap_genes_found_absolute_max_0100.eps'),
                       toolColors)

# Plot differences in whether the gene was found within the first 20 positions.
benchmarkQuintupleVenn(apply(positionResults <=20, 2, which),
                       sum(apply(positionResults > 20, 1, all, na.rm=T)),
                       paste0(imgExportDir, 'benchmarking_overlap_genes_found_absolute_max_0020.eps'),
                       toolColors)

# Plot differences in whether the gene was found within the first 20 positions.
benchmarkQuintupleVenn(apply(relativePositionResults <=0.2, 2, which),
                       sum(apply(relativePositionResults > 0.2, 1, all, na.rm=T)),
                       paste0(imgExportDir, 'benchmarking_overlap_genes_found_relative_max_0dot1.eps'),
                       toolColors)

postscript(paste0(imgExportDir, 'benchmarking_first_rank_phenotype_frequencies.eps'), width=10, height=6)
par(mar=c(11,4,4,0))
#barplot(t(phenotypeCountsWhenToolRankedFirst[apply(phenotypeCountsWhenToolRankedFirst, 1, function(x) { any(x>2) }),]),
#        col=colors, las=2, cex.names=0.5, space=0,
#        main="input phenotype frequencies when a tool ranked first in finding the gene\n(frequency > 2 for a single tool)",
#        ylab="phenotype input frequency")
dataToPlot <- t(phenotypeCountsWhenToolRankedFirst[which(phenotypeTotals > 4),])
barplot(dataToPlot,
        col=toolColors, las=2, cex.names=0.5, space=0,
        main="input phenotype frequencies when a tool ranked first in finding the gene\n(total frequency > 3)", # excludes input phenotypes with only NA
        ylab="phenotype input frequency")
legend(0,80, colnames(phenotypeCountsWhenToolRankedFirst),
       fill=toolColors)
dev.off()

postscript(paste0(imgExportDir, 'benchmarking_first_rank_phenotype_frequencies_relative.eps'), width=10, height=6)
par(mar=c(11,4,4,9))
barplot(prop.table(dataToPlot, 2),
        col=toolColors, las=2, cex.names=0.5, space=0, border="white",
        main="relative input phenotype frequencies when a tool ranked first in finding the gene\n(only phenotypes with total frequency > 3)", # excludes input phenotypes with only NA
        ylab="phenotype input frequency")
par(xpd=TRUE) # no clipping for drawing outside plot
legend(90,1, colnames(phenotypeCountsWhenToolRankedFirst),
       fill=toolColors)
dev.off()
