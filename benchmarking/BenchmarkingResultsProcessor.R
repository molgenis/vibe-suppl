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
# benchmarkDataRow - a list (row from benchmarkData given by
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

calculateTotalGenesFound <- function(benchmarkResults) {
  unlist(lapply(sapply(benchmarkResults[,"suggested_genes"], strsplit, split=","), length), use.names=FALSE)
}

sortRows <- function(benchmarkResults) {
  benchmarkResults[order(as.numeric(rownames(benchmarkResults))), , drop=FALSE]
}

########
# Name:
# drawBenchmarkQuintupleVenn
#
# Description:
# Draws a venn diagram using data as input and writes text outside of it indicating the number from outside to show
# how many were not included by any of the items from the venn diagram.
#
# Input:
# data - A list of vectors (e.g., integers, chars), with each component corresponding to a separate circle in the Venn
#        diagram (direct citation from venn.diagram() from the library VennDiagram).
# outside - a number indicating how many hits were outside the venn diagram.
#
# Output:
#
########
drawBenchmarkQuintupleVenn <- function(data, outside) {
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
vibe.20180503 <- readResultFile("results/vibe_2018-05-02.tsv")

# Sorts benchmark results so that row order is identical.
amelie <- sortRows(amelie)
geneNetwork <- sortRows(geneNetwork)
phenomizer <- sortRows(phenomizer)
phenotips <- sortRows(phenotips)
vibe.20180503 <- sortRows(vibe.20180503)



###
### Data calculations.
###

# Calculate absolute positions. This is done by processing benchmarkData row-by-row
# and therefore the LOVD row order of positionResults is equal to that of benchmarkData.
positionResults <- data.frame(amelie=resultsPositionCalculator(benchmarkData, amelie),
                              geneNetwork=resultsPositionCalculator(benchmarkData, geneNetwork),
                              phenomizer=resultsPositionCalculator(benchmarkData, phenomizer),
                              phenotips=resultsPositionCalculator(benchmarkData, phenotips),
                              vibe.20180503=resultsPositionCalculator(benchmarkData, vibe.20180503))

# Calculate total number of genes.
totalResults <- data.frame(amelie=calculateTotalGenesFound(amelie),
                           geneNetwork=calculateTotalGenesFound(geneNetwork),
                           phenomizer=calculateTotalGenesFound(phenomizer),
                           phenotips=calculateTotalGenesFound(phenotips),
                           vibe.20180503=calculateTotalGenesFound(vibe.20180503),
                           row.names=rownames(amelie))
totalResults <- totalResults[benchmarkData$lovd,]


# Calculate relative positions.
relativePositionResults <- positionResults / totalResults

# Count missing values.
resultsWithNa <- positionResults[apply(apply(positionResults, c(1,2), is.na), 1, any),]
naCounts <- apply(apply(resultsWithNa, c(1,2), is.na), 2, sum)

### Calculate absolute tool rankings.
toolRankingResults <- sapply(apply(positionResults, 1, sort, na.last=NA), names)
toolRankingResults <- sapply(1:5, function(x) { sapply(toolRankingResults, '[', x) })
colnames(toolRankingResults) <- c("first", "second", "third", "fourth", "fifth")

toolRankingResultCounts <- apply(toolRankingResults, 2, function(x) {
  table(factor(x, levels=colnames(positionResults)))
  })
toolRankingResultCounts <- cbind(toolRankingResultCounts, "no hit"=naCounts)

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
  sapply(rownames(toolRankingResultCounts), function(toolName) {
    table(unlist(
      sapply(benchmarkData[which(toolRankingResults[,"first"]==toolName),5],
             strsplit, ";"),
      use.names=FALSE))
    })

# Adds the results from foundCounts to phenotypeCountsWhenToolRankedFirst.
for(toolName in colnames(phenotypeCountsWhenToolRankedFirst)) {
  phenotypeCountsWhenToolRankedFirst[names(foundCounts[[toolName]]),toolName] <- foundCounts[[toolName]]
}
rm(foundCounts)



###
### Data info.
###

# Check whether the rows of all benchmarks are identically ordered.
rownames(amelie) == rownames(geneNetwork) &&
  rownames(geneNetwork) == rownames(phenomizer) &&
  rownames(phenomizer) == rownames(phenotips) &&
  rownames(phenotips) == rownames(vibe.20180503)

# Check whether there is any phenotype set for which no single gene was found.
any(is.na(totalResults))

# Check whether tool ranking return values indicating they are complete (and equal!).
apply(toolRankingResultCounts, 1, sum)



###
### Plotting figures.
###

# Boxplot comparing absolute positions of genes.
postscript(paste0(imgExportDir, 'benchmarking_gene_position_absolute.eps'), width=7, height=8)
yAxisMax <- ceiling(max(positionResults, na.rm=T)/100)*100
boxplot(positionResults, xaxt='n',yaxt='n', pch="")
abline(h=1, col="gray92")
abline(h=seq(100, yAxisMax, 100), col="gray92")
boxplot(positionResults, las=1, pch=20, yaxt='n',
        main="position of relevant genes among different tools",
        xlab="tool", ylab="position", col="white", add=T)
axis(2, las=1, at=1)
axis(2, las=1, at=seq(500, yAxisMax,500))
dev.off()

# Boxplot comparing number of suggested genes found.
postscript(paste0(imgExportDir, 'benchmarking_gene_position_total_hits.eps'), width=7, height=8)
yAxisMax <- ceiling(max(totalResults, na.rm=T)/500)*500
boxplot(totalResults, xaxt='n',yaxt='n', pch="")
abline(h=1, col="gray92")
abline(h=seq(500, yAxisMax, 500), col="gray92")
boxplot(totalResults, las=1, pch=20, yaxt='n',
        main="the number of output genes among different tools",
        xlab="tool", ylab="number of genes", col="white", add=T)
axis(2, las=1, at=1)
axis(2, las=1, at=seq(2000, yAxisMax,2000))
dev.off()

# Boxplot comparing relative positions of genes (absolute position / number of suggested genes found).
postscript(paste0(imgExportDir, 'benchmarking_gene_position_relative.eps'), width=7, height=8)
boxplot(relativePositionResults, xaxt='n',yaxt='n', pch="")
abline(h=seq(0, 1, 0.05), col="gray92")
boxplot(relativePositionResults, las=1, pch=20,
        main="relative position of relevant genes among different tools",
        xlab="tool", ylab="relative position", col="white", add=T)
dev.off()

# Barplot showing the tool rankings.
postscript(paste0(imgExportDir, 'benchmarking_tool_ranking.eps'), width=7, height=4)
colors <- rev(brewer.pal(6, 'RdYlGn'))
barplot(t(toolRankingResultCounts), col=colors, las=1,
        main="gene position ranked among the tools")
par(xpd=TRUE) # no clipping for drawing outside plot
legend(0,-80, colnames(toolRankingResultCounts),
       fill=colors, ncol=ncol(toolRankingResultCounts))
dev.off()

# Venn diagram basic settings.
colors <- brewer.pal(ncol(positionResults), 'Set3')

# Plot differences in whether the gene was found.
cairo_ps(paste0(imgExportDir, 'benchmarking_overlap_genes_found_absolute.eps'), width=7, height=7)
drawBenchmarkQuintupleVenn(apply(!is.na(positionResults), 2, which),
                           sum(apply(is.na(positionResults), 1, all)))
dev.off()

# Plot differences in whether the gene was found within the first 100 positions.
cairo_ps(paste0(imgExportDir, 'benchmarking_overlap_genes_found_absolute_max_100.eps'), width=7, height=7)
drawBenchmarkQuintupleVenn(apply(positionResults <=100, 2, which),
                           sum(apply(positionResults > 100, 1, all, na.rm=T)))
dev.off()

# Plot differences in whether the gene was found within the first 20 positions.
cairo_ps(paste0(imgExportDir, 'benchmarking_overlap_genes_found_absolute_max_020.eps'), width=7, height=7)
drawBenchmarkQuintupleVenn(apply(positionResults <=20, 2, which),
                           sum(apply(positionResults > 20, 1, all, na.rm=T)))
dev.off()

# Plot differences in whether the gene was found within the first 20 positions.
cairo_ps(paste0(imgExportDir, 'benchmarking_overlap_genes_found_relative_max_0dot1.eps'), width=7, height=7)
drawBenchmarkQuintupleVenn(apply(relativePositionResults <=0.2, 2, which),
                           sum(apply(relativePositionResults > 0.2, 1, all, na.rm=T)))
dev.off()

postscript(paste0(imgExportDir, 'benchmarking_first_rank_phenotype_frequencies.eps'), width=10, height=6)
par(mar=c(10,4,4,2))
barplot(t(phenotypeCountsWhenToolRankedFirst[apply(phenotypeCountsWhenToolRankedFirst, 1, function(x) { any(x>2) }),]),
        col=colors, las=2, cex.names=0.5, space=0,
        main="input phenotype frequencies when a tool ranked first in finding the gene\n(frequency > 2 for a single tool)",
        ylab="phenotype input frequency")
legend(0,80, colnames(phenotypeCountsWhenToolRankedFirst),
       fill=colors)
dev.off()
par(oldPar)
