########
# Name:
# BenchmarkResultsGenerics.R
#
# Description:
# Generic functions, data and such for the different benchmarking result scripts.
# 
# Important:
# Be sure to adjust the items in the config section before running the script!
########

##################
### Config     ###
##################

# Directory storing script input:
# - benchmark_input/benchmark_data.tsv
# - results/amelie.tsv
# - results/gene_network.tsv
# - results/phenomizer.tsv
# - results/phenotips.tsv
# - results/vibe_2018-07-06_gda_max.tsv
# - results/vibe_2018-07-06_dsi.tsv
# - results/vibe_2018-07-06_dpi.tsv
baseDir <- '~/vibe_benchmarking_r_input/'

# Directory to write the output images to.
# Be sure this directory contains the following subfolders:
# - tools
# - vibe_sorting_algorithms
# - vibe_adjustments
baseImgExportDir <- '~/vibe_benchmarking_r_output/'





##################
### Libraries  ###
##################

library(VennDiagram)
library(RColorBrewer)
library(dplyr)
#library(scales) # imported through ggplot2
library(reshape2)
library(ggplot2)
library(ggdendro)
library(viridis)
#library(heatmaply) # code commented due to issues
library(cowplot)





##################
### Functions  ###
##################

########
# Name:
# initializeGraphicsDevice
#
# Description:
# A simple wrapper for creating a graphicsdevice that prepends the directory
# and appends the file extension.
#
# Input:
# fileName - Filename (excluding file extension) to be used for storage.
# width	- The width of the device in inches.
# height - The height of the device in inches.
#
# Output:
#
# Important:
# To allow for easy adjustments between different output engines, inches is
# used as default output width/height and for engines using pixels a simple
# multiplier is used.
########
initializeGraphicsDevice <- function(fileName, width, height) {
  #cairo_ps(paste0(imgExportDir, fileName, ".eps"), width=width, height=height)
  pdf(paste0(imgExportDir, fileName, ".pdf"), width=width, height=height)
}

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
# readTimesFile
#
# Description:
# Reads in a single file containing benchmarking times.
#
# Input:
# filePath -  The path to the file to be loaded.
#
# Output:
# A table from the result data. Though as the input is expected to only contain
# 2 columns and one of these is used as row names, type becomes a list.
########
readTimesFile <- function(filePath) {
  read.table(paste0(baseDir, filePath),
             header=T, sep="\t", colClasses=c("character","double"),
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
  return(benchmarkResults[order(as.numeric(rownames(benchmarkResults))), , drop=FALSE])
}

########
# Name:
# calculateRankings
#
# Description:
# Converts a dataframe with tool columns and cells as the benchmark score to
# a matrix that contains columns for the ranking and the cell values with the
# toolname that has that ranking. For example: If input contains 2 columns
# (X & Y) where column Y has the lowest value for a row, the output will contain
# Y as cell value for the column named "first".
#
# IMPORTANT: Currently supports up to 5 columns as input (up to "fifth").
#
# Input:
# data - A dataframe which needs to be sorted. The columns should contain the
#        tools that need to be ranked. The rows the individual benchmark tests.
#        Each cell then contains the benchmark position for the given tool for
#        that benchmark test.
#
# Output:
# A character matrix with as columns the positions (first, second, etc.), as
# rows the individual tests and as cell values the tool names.
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
# Converts the output from calculateRankings() to a integer matrix containing
# how often each tool scored which position. Furthermore, adds the NA counts
# to this matrix as well.
#
# IMPORTANT: Make sure that the output from rankingResults() used the same data
#            as input!!!
#
# Input:
# data - A dataframe which needs to be sorted. The columns should contain the
#        tools that need to be ranked. The rows the individual benchmark tests.
#        Each cell then contains the benchmark position for the given tool for
#        that benchmark test.
# rankingResults - The output from calculateRankings() when using THE SAME data
#                  as input!
#
# Output:
# An integer matrix with as rows the different tools and as columns the
# rankings ("first", "second", "third" and so forth with "no hit" for if a tool
# did not find the match at all).
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
# Generates a boxplot using the data as input.
#
# Input:
# data - The dataframe to create a plot from.
# fileName - Filename (excluding file extension) to be used for storage.
# firstYValue - The y-axis starting position.
# yAxisFreq - The frequency for the y-axis labels.
# yAxisAblineFreq - The frequency of the vertical lines in the plot.
# main - The plot title.
# ylab - The y-axis label.
# parMar - the margin of the graphical parameters. DEFAULT=c(4.5,4.5,0.5,0.5)
# ylabMgp - Adjusts the mgp setting for ylab (see ?par). Only the first parameter
#           is relevant. DEFAULT=c(3,1,0)
#
# Output:
#
########
toolValuesBoxPlot <- function(data, fileName, firstYValue, yAxisFreq,
                              yAxisAblineFreq, main, ylab,
                              parMar=c(4.5,4.5,0.5,0.5), ylabMgp=c(3,1,0)) {
  initializeGraphicsDevice(fileName, width=7, height=8)
  par(mar=parMar)
  
  yAxisMax <- ceiling(max(data, na.rm=T)/yAxisFreq)*yAxisFreq
  boxplot(data, ylim=c(firstYValue,yAxisMax), xaxt='n',yaxt='n', pch="")
  abline(h=seq(firstYValue, yAxisMax+firstYValue, yAxisAblineFreq), col="gray85")
  boxplot(data, las=1, pch=20, yaxt='n',
          main=main, xlab="tool", ylab="", col="white", add=T)
  title(ylab=ylab, mgp=ylabMgp)
  axis(2, las=1, at=seq(firstYValue, yAxisMax+firstYValue, yAxisFreq),
       labels=format(seq(firstYValue, yAxisMax+firstYValue, yAxisFreq),
                     big.mark=","))
  
  dev.off()
}

########
# Name:
# toolRankingBarplot
#
# Description:
# Generates a barplot based on the output from calculateRankingCounts().
#
# Input:
# data - The data to be plotted (in general: the object generated from
#        calculateRankingCounts() ).
# fileName - Filename (excluding file extension) to be used for storage.
# main - The plot title.
# ylab - The y-axis label.
# legendPos - A vector with 2 positions used for placing the legend().
#             The first item from the vector will be used for the x-axis
#             positioning. The second item from the vector will be used for the
#             y-axis positioning.
# parMar - the margin of the graphical parameters. DEFAULT=c(5.0,4.5,1.0,0.5)
#
# Output:
#
########
toolRankingBarplot <- function(data, fileName, main, ylab, legendPos, parMar=c(5.0,4.5,1.0,0.5)) {
  GreenToRedColors <- rev(brewer.pal(nrow(data), 'RdYlGn'))
  
  yAxisFactor <- 20
  totalValues <- max(apply(data, 2, sum))
  yAxisMax <- ceiling(totalValues/yAxisFactor)*yAxisFactor
  
  initializeGraphicsDevice(fileName, width=7, height=8)
  par(mar=parMar)
  
  barplot(data, col=GreenToRedColors, las=1, yaxt="n",
          main=main, ylab=ylab, ylim=c(0,yAxisMax))
  axis(2, at=seq(0, yAxisMax, yAxisFactor), las=1)
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
# fileName - Filename (excluding file extension) to be used for storage.
# colors - The colors to be used for the venn diagram.
#
# Output:
#
########
benchmarkQuintupleVenn <- function(data, outside, fileName, colors) {
  initializeGraphicsDevice(fileName, width=7, height=7)
  
  grid.draw(venn.diagram(data, NULL,
                         fontfamily="Helvetica", main.fontfamily="Helvetica",
                         sub.fontfamily="Helvetica", cat.fontfamily="Helvetica",
                         fill=colors, cat.just=list(c(0.5,1),
                                                    c(-0.5,-4),
                                                    c(0,0),
                                                    c(1,1),
                                                    c(3.5,-5))
  ))
  grid.text(paste("none\n", outside), 0.1, 0.2)
  
  dev.off()
}

########
# Name:
# hitPositionsScatterplot
#
# Description:
# Creates a scatterplot based on the data. A vertical line (Y=X) is also included.
#
# Input:
# data - The data to be plott. Should contain 3 columns: the first column should
#        contain the x-axis values, the second column the y-axis values and the
#        third column should contain the different groups (are assigned
#        different colors).
# fileName - Filename (excluding file extension) to be used for storage.
# colors - The colors to be used within the plot (each group gets a different
#          color)
# xAxisFreq - The frequency for the x-axis labels.
# xAxisAblineFreq - The frequency of x-axis lines.
# yAxisFreq - The frequency for the y-axis labels.
# yAxisAblineFreq - The frequency of y-axis lines.
# xlab - The x-axis label.
# ylab - The y-axis label.
# parMar - the margin of the graphical parameters. DEFAULT=c(4.5,4.5,0.5,0.5)
# ylabMgp - Adjusts the mgp setting for ylab (see ?par). Only the first parameter
#           is relevant. DEFAULT=c(3,1,0)
#
# Output:
#
########
hitPositionsScatterplot <- function(data, fileName, colors, xAxisFreq,
                                    xAxisAblineFreq, yAxisFreq, yAxisAblineFreq,
                                    xlab, ylab, parMar=c(4.5,4.5,0.5,0.5),
                                    ylabMgp=c(3,1,0)) {
  # Renames colnames for universal usage.
  colnames(data) <- c("x", "y", "group")
  
  # Replaces Inf with NA as these values cannot be plotted (but could cause
  # issues if left in the data).
  data[which(abs(data$x) == Inf), "x"] <- NA
  data[which(abs(data$y) == Inf), "y"] <- NA
  
  # Removes values that have an NA as either x or y coordinate.
  data <- data[apply(data[,c("x", "y")], 1, function(x) {!any(is.na(x))}),]
  
  # Adjustments for if the group column has more levels than mentioned by the
  # actual rows.
  colors <- colors[c(unique(data[,3]))]
  data[,3] <- droplevels(data[,3])
  
  # Calculates highest x & y value.
  xAxisMax <- ceiling(max(data[,1])/xAxisAblineFreq)*xAxisAblineFreq
  yAxisMax <- ceiling(max(data[,2])/yAxisAblineFreq)*yAxisAblineFreq
  
  # Recurring usage.
  toolNames <- levels(data[,3])
  
  # Plot total genes against found position.
  initializeGraphicsDevice(fileName, width=10, height=6)
  par(mar=parMar)
  
  plot(1, las=1, type="n", bty="u",
       xlab=xlab,
       ylab="",
       xaxt="n",
       yaxt="n",
       xlim=c(0, xAxisMax),
       ylim=c(0, yAxisMax))
  axis(1, las=1, at=seq(0, xAxisMax, xAxisFreq),
       labels=format(seq(0, xAxisMax, xAxisFreq), big.mark=","))
  axis(2, las=1, at=seq(0, yAxisMax, yAxisFreq),
       labels=format(seq(0, yAxisMax, yAxisFreq), big.mark=","))
  title(ylab=ylab, mgp=ylabMgp)
  abline(h=seq(0,yAxisMax, yAxisAblineFreq),
         v=seq(0,xAxisMax, xAxisAblineFreq), col="gray85")
  abline(0,1, col="gray85") # indicates where position limit is (position <= hits)
  legend("topleft", c(toolNames, "combined", "patient case",
                      "mean", "linear model"),
         col=c(colors, "red", rep("black", 3)),
         pch=c(rep(15, length(toolNames)+1), 20, 25, NA),
         lty=c(rep(NA, length(toolNames)+3), 1),
         bg="white")
  points(y~x, data, col=alpha(colors[group], 0.4), pch=20)
  
  # Calculates data per group.
  groupData <- data %>% group_by(group) %>% summarise(x.mean = mean(x),
                                                      y.mean = mean(y),
                                                      x.min = min(x),
                                                      x.max = max(x),
                                                      y.min = min(y),
                                                      y.max = max(y))
  
  # Linear model of all data.
  clip(min(data$x), max(data$x), min(data$y), max(data$y))
  abline(lm(y~x, data), col="red", lwd=2)
  
  # Linear model per tool.
  groupLinearModels <- data %>% group_by(group) %>% do(model = lm(y~x, data = .))
  for(i in 1:nrow(groupLinearModels)) {
    # Clips range for abline to the limits of the specific tool.
    clip(groupData$x.min[[i]],
         groupData$x.max[[i]],
         groupData$y.min[[i]],
         groupData$y.max[[i]])
    # A try is required for cases where a tool always returns the exact number
    # of output genes (so x values are always the same).
    try(abline(groupLinearModels$model[[i]], col=colors[i], lwd=2))
  }
  
  # Resets clip region.
  do.call("clip", as.list(par("usr")))
  
  # Adds mean per group to plot.
  points(y.mean~x.mean, groupData, bg=colors[group], pch=25, col="black")
  
  dev.off()
}

########
# Name:
# plotToolComaprison
#
# Description:
# Plots 2 tools against each other on how they compare for each row from the
# input dataframe.
#
# Input:
# data - The data to be plotted. Should contain 2 columns containing the different
#        tools that are compared.
# fileName - Filename (excluding file extension) to be used for storage.
# xyMax - The plot max of both the x-axis and y-axis.
# xlab - The x-axis label.
# ylab - The y-axis label.
# xyMin - The plot min of both the x-axis and y-axis. DEFAULT=0
# naHitsAdjust - Axis adjustment for plotting values outside the plot
#                (indicating only 1 of the tools has a result while the other
#                has an NA there). DEFAULT=0.3
# parMar - the margin of the graphical parameters. DEFAULT=c(4.5,4.5,1.5,1.5)
#
# Output:
#
########
plotToolComaprison <- function(data, fileName, xyMax, xlab, ylab, xyMin=0,
                               naHitsAdjust=0.3, parMar=c(4.5,4.5,1.5,1.5)) {
  col1ValuesWithNaCol2 <- data[is.na(data[,2]),1]
  col2ValuesWithNaCol1 <- data[is.na(data[,1]),2]
  
  initializeGraphicsDevice(fileName, width=6, height=6)
  par(mar=parMar)
  plot(data[,1], data[,2], las=1, pch=20,
       xlim=c(xyMin,xyMax), ylim=c(xyMin,xyMax), xlab=xlab, ylab=ylab,
       # Makes sure grid & legend are drawn before the data.
       panel.first=c(grid(col="gray85", lty=1),
                     legend("topleft",
                            c("patient case", "mean", "Y=X", "linear model"),
                            pch=c(20, 20, NA, NA),
                            lty=c(NA, NA, 1, 1),
                            col=c("black", "orange", "red", "orange"),
                            bg="white")))
  abline(0,1, col="red")
  abline(lm(data[,2]~data[,1]), col="orange")
  par(xpd=TRUE) # no clipping for drawing outside plot
  points(col1ValuesWithNaCol2,
         rep(xyMax + naHitsAdjust, length(col1ValuesWithNaCol2)), pch=20)
  points(rep(xyMax + naHitsAdjust, length(col2ValuesWithNaCol1)),
         col2ValuesWithNaCol1, pch=20)
  points(mean(data[,1], na.rm=TRUE), 
         mean(data[,2], na.rm=TRUE), col="orange", pch=20)
  dev.off()
}

########
# Name:
# plotMatchesFoundWihinRangeCutoff
#
# Description:
# Plots how many hits were found within a defined cutoff. The data should
# contain the number/fraction of found hits within the cutoff defined by the
# row names.
#
# Input:
# data - The data to be plotted. Should contain a column per tool to plot and
#        row names indicating the x-axis values.
# fileName - Filename (excluding file extension) to be used for storage.
# xlab - The x-axis label.
# ylab - The y-axis label.
# colors - The colors to be used within the plot (each group gets a different
#          color)
# xLabelsStepSize - The step size of the x-axis labels.
# addAbline - Whether abline(0,1) should be added to plot. DEFAULT=TRUE
# parMar - the margin of the graphical parameters. DEFAULT=c(4.5,4.5,0.5,0.5)
#
# Output:
#
########
plotMatchesFoundWihinRangeCutoff <- function(data, fileName, xlab, ylab, colors,
                                             xLabelsStepSize, addAbline=TRUE,
                                             parMar=c(4.5,4.5,0.5,0.5)) {
  highestUsedCutoffWithinData <- as.numeric(rownames(data)[nrow(data)])
  
  initializeGraphicsDevice(fileName, width=10, height=5)
  par(mar=parMar)
  # Generates empty plot to use.
  plot(1, type="n", las=1, xlim=c(0, highestUsedCutoffWithinData), ylim=c(0,1),
       xlab=xlab, ylab=ylab, xaxt="n", yaxt="n")
  axis(1, at=seq(0,highestUsedCutoffWithinData,xLabelsStepSize), las=1)
  axis(2, at=seq(0,1,0.1), las=1)
  
  # Adds abline if requested.
  if(addAbline) {abline(0,1, col="gray50")}
  
  # Adds data lines.
  for(x in 1:ncol(data)) {
    lines(rownames(data), data[,x], col=colors[x])
  }
  
  # Adds legend.
  if(addAbline) {
    legend("bottomright", c(colnames(data), "Y=X"),
           pch=c(rep(15, length(colnames(data))), NA),
           lty=c(rep(NA, length(colnames(data))), 1),
           col=c(colors, "gray50"))
  } else {
    legend("bottomright", colnames(data), pch=15, col=colors)
  }
  
  dev.off()
}

########
# Name:
# ggplotLegend
#
# Description:
# Retrieves legend from a ggplot.
# Based on https://stackoverflow.com/a/12041779
#
# Input:
# ggplot - The plot to retrieve the legend from.
#
# Output:
# The legend from the plot.
########
ggplotLegend <- function(ggplot){ 
  tmp <- ggplot_gtable(ggplot_build(ggplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
}

########
# Name:
# toolScoresHistogramWithDensity
#
# Description:
# Plots histograms
#
# Input:
# data  - The data to be plotted. Should one-dimensional (no data frame).
# fileName - Filename (excluding file extension) to be used for storage.
# xlab  - The x-axis label.
# xlim  - vector with two values describing the lower and upper bound for the
#         x-axis.
# ylimMax - The max y-axis value.
#
# Output:
# 
########
toolScoresHistogram <- function(data, fileName, xlab, xlim, ylimMax) {
  initializeGraphicsDevice(fileName, width=10, height=5)
  par(mar=c(5.5,5.5,0.5,0.5), mgp=c(4,1,0))
  hist(data, main="", xlab=xlab, ylab="probability density", freq=F, las=1, breaks=20,
       xlim=xlim, ylim=c(0,ylimMax), col="steelblue4")
  lines(density(data, na.rm=T), col="red", lwd=2)
  legend(xlim[1], ylimMax, "density curve", lty=1, col="red")
  dev.off()
}


########
# Name:
# calculateTotalGeneFrequencies
#
# Description:
# Calculates the total frequencies of genes among the whole benhmark output.
#
# Input:
# benchmarkResults  - Results from a benchmark.
#
# Output:
# A table containing the genes with their frequency counts.
########
calculateTotalGeneFrequencies <- function(benchmarkResults) {
  convertedData <- sapply(benchmarkResults$suggested_genes, strsplit, split=",", USE.NAMES = FALSE)
  convertedData <- do.call(c, convertedData)
  geneCounts <- sort(table(convertedData), decreasing = TRUE)
  return(geneCounts)
}

########
# Name:
# calculateLovdCountsForGenesWithMultipleOccurencesInSingleLovd
#
# Description:
# Calculates for genes that occur multiple times in the output for a singe LOVD,
# for how many LOVDs this is the case.
#
# Input:
# benchmarkResults  - Results from a benchmark.
#
# Output:
# A vector that contains the genes as names and as values for how many LOVDs this gene
# was found multiple times within the output.
########
calculateLovdCountsForGenesWithMultipleOccurencesInSingleLovd <- function(benchmarkResults) {
  genesPerLovd <- sapply(benchmarkResults$suggested_genes, function(suggestedGenes) {
    geneFrequencies <- table(strsplit(suggestedGenes, ","))
    return(names(which(geneFrequencies > 1)))
  })
  numberOfLovdsGenesHaveMultipleOccurencesIn <- table(unlist(genesPerLovd))
  return(numberOfLovdsGenesHaveMultipleOccurencesIn)
}

########
# Name:
# plotTotalGeneFrequencyfunction
#
# Description:
# 
#
# Input:
# 
#
# Output:
# A vector that contains the genes as names and as values for how many LOVDs this gene
# was found multiple times within the output.
########
plotTotalGeneFrequencyfunction <- function(dataToPlot, totalCases, yLimMax, title) {
  barplot(dataToPlot, las=2, col="steelblue4", space=0, border=NA, ylim=c(0,yLimMax),
          names.arg=NA, main=title,
          xlab="genes ordered by frequency\n(too many for individual names, count shows number of genes)", ylab='number of patient cases gene is found in output')
  abline(h=totalCases, col="red")
  legend(length(dataToPlot)*0.5,yLimMax, c('number of LOVDs'), lty=1, col="red")
  axis(1, at=c(length(dataToPlot)))
}

########
# Name:
# plotGenesWithMultipleOccurencesInSingleLovds
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
plotGenesWithMultipleOccurencesInSingleLovds <- function(dataToPlot, yLimMax, title) {
  xLab <- 'genes that occured multiple times in a single LOVD output'
  yLab <- 'n cases gene occured multiple times within a single LOVD output'
  
  if(length(dataToPlot) == 0) {
    barplot(0, las=2, ylim=c(0,yLimMax), main=title, ylab=yLab, col=NA, border=NA)
  }
  else {
    barplot(dataToPlot, las=2, col="steelblue4", ylim=c(0,yLimMax), main=title, ylab=yLab)
  }
  title(xlab = xLab, line = 6)
}




##################
###    Code    ###
##################

# Defaults
oldPar <- par()
setEPS() # Sets EPS engine for writing images.

# Load benchmark data.
benchmarkData <- read.table(paste0(baseDir,"benchmark_data.tsv"), header=T,
                            sep="\t",colClasses=c(rep("character", 3),
                                                  "factor", "character"))
