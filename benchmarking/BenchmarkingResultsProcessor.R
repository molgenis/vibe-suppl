########
# Name:
# BenchmarkingResultsProcessor.R
#
# Description:
# Generates plots from the result data.
# 
# Important:
# Be sure to adjust the items in the config section before running the script!
########

##################
### Config     ###
##################

# Directory storing script input:
# - benchmark_data.tsv
# - results/amelie.tsv
# - results/gene_network.tsv
# - results/phenomizer.tsv
# - results/phenotips.tsv
# - results/vibe_2018-07-06_gda_max.tsv
# - results/vibe_2018-07-06_dsi.tsv
# - results/vibe_2018-07-06_dpi.tsv
baseDir <- '~/programming/projects/vibe/data/benchmarking/'

# Directory to write the output images to.
imgExportDir <- '~/Documents/afstuderen/verslag/img/'





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
# A simple wrapper for cairo_ps() that prepends the directory and appends the
# file extension.
#
# Input:
# fileName - Filename (excluding file extension) to be used for storage.
# width	- The width of the device in inches. -> see ?cairo_ps()
# height - The height of the device in inches. -> see ?cairo_ps()
#
# Output:
#
########
initializeGraphicsDevice <- function(fileName, width, height) {
  cairo_ps(paste0(imgExportDir, fileName, ".eps"), width=width, height=height)
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
  if(addAbline) {abline(0,1, col="gray70")}
  
  # Adds data lines.
  for(x in 1:ncol(data)) {
    lines(rownames(data), data[,x], col=colors[x])
  }
  
  # Adds legend.
  if(addAbline) {
    legend("bottomright", c(colnames(data), "Y=X"),
           pch=c(rep(15, length(colnames(data))), NA),
           lty=c(rep(NA, length(colnames(data))), 1),
           col=c(colors, "gray70"))
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





##################
###    Code    ###
##################

###
### Basic settings & data basic preparation.
###

# Defaults
oldPar <- par()
setEPS() # Sets EPS engine for writing images.
toolColors <- brewer.pal(5, 'Set2') # Adjust if changing number of tools!!!

# Load data.
benchmarkData <- read.table(paste0(baseDir,"benchmark_data.tsv"), header=T,
                            sep="\t",colClasses=c(rep("character", 3),
                                                  "factor", "character"))

amelie <- readResultFile("results/amelie.tsv")
gado <- readResultFile("results/gene_network.tsv")
phenomizer <- readResultFile("results/phenomizer.tsv")
phenotips <- readResultFile("results/phenotips.tsv")

vibe.gda_max <- readResultFile("results/vibe_2018-07-06_gda_max.tsv")
vibe.dsi <- readResultFile("results/vibe_2018-07-06_dsi.tsv")
vibe.dpi <- readResultFile("results/vibe_2018-07-06_dpi.tsv")

# Sorts benchmark results so that row order is identical.
amelie <- sortRows(amelie)
gado <- sortRows(gado)
phenomizer <- sortRows(phenomizer)
phenotips <- sortRows(phenotips)

vibe.gda_max <- sortRows(vibe.gda_max)
vibe.dsi <- sortRows(vibe.dsi)
vibe.dpi <- sortRows(vibe.dpi)

# Loads in best score per tool.
toolGeneScores <- read.table(paste0(baseDir, "tool_gene_scores.tsv"),
                             header=T, sep="\t")
toolGeneScores <- toolGeneScores[,3:5]



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
### Calculations for comparison of different tools (with best scoring vibe result).
###

# Sets which vibe algorithm should be used for further processing.
vibe <- vibe.gda_max


# Calculate absolute positions. This is done by processing benchmarkData row-by-row
# and therefore the LOVD row order of positionResults is equal to that of benchmarkData.
positionResults <- data.frame(amelie=resultsPositionCalculator(benchmarkData, amelie),
                              gado=resultsPositionCalculator(benchmarkData, gado),
                              phenomizer=resultsPositionCalculator(benchmarkData, phenomizer),
                              phenotips=resultsPositionCalculator(benchmarkData, phenotips),
                              vibe=resultsPositionCalculator(benchmarkData, vibe))

# Calculate total number of genes.
totalResults <- data.frame(amelie=calculateTotalGenesFound(amelie),
                           gado=calculateTotalGenesFound(gado),
                           phenomizer=calculateTotalGenesFound(phenomizer),
                           phenotips=calculateTotalGenesFound(phenotips),
                           vibe=calculateTotalGenesFound(vibe),
                           row.names=rownames(amelie))

# Replicates some of the totalResults so that size is equal to positionResults.
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
                                                 gado=rep(0,length(allPhenotypeNames)),
                                                 phenomizer=rep(0,length(allPhenotypeNames)),
                                                 phenotips=rep(0,length(allPhenotypeNames)),
                                                 vibe=rep(0,length(allPhenotypeNames)),
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

# Calculates how many of the genes are found when only looking within specific
# absolute cutoffs of the total number of hits available.
highestPossibleResultForTool <- apply(totalResults, 2, max)
cutoffRanges <- seq(1, max(highestPossibleResultForTool), 1) # how specific looking for cutoffs
genesFoundWithinAbsoluteCutoff <- apply(positionResults, 2, function(toolAbsoluteScores, range=cutoffRanges) {
  sapply(range, function(cutoff, scores=toolAbsoluteScores) {
    length(which(scores < cutoff))
  })
})
rownames(genesFoundWithinAbsoluteCutoff) <- cutoffRanges

# Replaces values with NA for positions if higher than highest available
# totalResult for that tool.
for(x in names(highestPossibleResultForTool)) {
  if(highestPossibleResultForTool[x] < nrow(genesFoundWithinAbsoluteCutoff)) {
    genesFoundWithinAbsoluteCutoff[(highestPossibleResultForTool[x]+1):nrow(genesFoundWithinAbsoluteCutoff),x] <- NA
  }
}

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
### Data info.
###

# Check whether the rows of all benchmarks are identically ordered.
rownames(amelie) == rownames(gado) &&
  rownames(gado) == rownames(phenomizer) &&
  rownames(phenomizer) == rownames(phenotips) &&
  rownames(phenotips) == rownames(vibe)

# Check whether there is any phenotype set for which no single gene was found.
any(is.na(totalResults))

# Check whether tool ranking return values indicating they are complete (and equal!).
apply(toolRankingCounts, 1, sum)

# Retrieves rows for which have an NA (except if only phenomizer has an NA).
naExceptPhenomizerOnly <- positionResults[apply(positionResults, 1,
                                                function(x) {
                                                  any(is.na(x[c(1,2,4,5)]))
                                                }),]

# Combines previous data with the lovd and gene information for each row.
naExceptPhenomizerOnly <- cbind(benchmarkData[rownames(naExceptPhenomizerOnly),
                                              c("lovd", "gene")],
                                naExceptPhenomizerOnly)



###
### Plotting vibe-only figures.
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
### Basic boxplots/barplots comparing the different tools.
###

# Boxplot comparing absolute positions of genes.
toolValuesBoxPlot(log10(positionResults),
                  'benchmarking_gene_position_absolute',
                  0, 1, 0.5,
                  "",
                  "absolute position of validation genes within the output genes (log10)")

# Boxplot comparing number of suggested genes found.
toolValuesBoxPlot(log10(totalResults),
                  'benchmarking_gene_position_total_hits',
                  0, 1, 0.5,
                  "",
                  "number of output genes (log10)")

# Boxplot comparing relative positions of genes (absolute position / number of suggested genes found).
toolValuesBoxPlot(relativePositionResults,
                  'benchmarking_gene_position_relative',
                  0.0, 0.2, 0.1,
                  "",
                  "relative position of validation genes within the output genes")

# Barplot showing the tool rankings.
toolRankingBarplot(t(toolRankingCounts),
                   'benchmarking_tool_ranking',
                   "", "patient cases", c(0,-23))



###
### Venn diagrams at different cutoffs comparing overlap in found genes between tools.
###

# Plot differences in whether the gene was found.
benchmarkQuintupleVenn(apply(!is.na(positionResults), 2, which),
                       sum(apply(is.na(positionResults), 1, all)),
                       'benchmarking_overlap_genes_found_absolute',
                       toolColors)

# Plot differences in whether the gene was found within the first 10000 positions.
benchmarkQuintupleVenn(apply(positionResults <=10000, 2, which),
                       sum(apply(positionResults > 10000, 1, all, na.rm=T)),
                       'benchmarking_overlap_genes_found_absolute_max_10000',
                       toolColors)

# Plot differences in whether the gene was found within the first 1000 positions.
benchmarkQuintupleVenn(apply(positionResults <=1000, 2, which),
                       sum(apply(positionResults > 1000, 1, all, na.rm=T)),
                       'benchmarking_overlap_genes_found_absolute_max_1000',
                       toolColors)

# Plot differences in whether the gene was found within the first 100 positions.
benchmarkQuintupleVenn(apply(positionResults <=100, 2, which),
                       sum(apply(positionResults > 100, 1, all, na.rm=T)),
                       'benchmarking_overlap_genes_found_absolute_max_100',
                       toolColors)

# Plot differences in whether the gene was found within the first 20 positions.
benchmarkQuintupleVenn(apply(positionResults <=20, 2, which),
                       sum(apply(positionResults > 20, 1, all, na.rm=T)),
                       'benchmarking_overlap_genes_found_absolute_max_20',
                       toolColors)

# Plot differences in whether the gene was found within the first 10 positions.
benchmarkQuintupleVenn(apply(positionResults <=10, 2, which),
                       sum(apply(positionResults > 10, 1, all, na.rm=T)),
                       'benchmarking_overlap_genes_found_absolute_max_10',
                       toolColors)

# Plot differences in whether the gene was found within the first 40% of output genes.
benchmarkQuintupleVenn(apply(relativePositionResults <=0.4, 2, which),
                       sum(apply(relativePositionResults > 0.4, 1, all, na.rm=T)),
                       'benchmarking_overlap_genes_found_relative_max_0dot4',
                       toolColors)

# Plot differences in whether the gene was found within the first 30% of output genes.
benchmarkQuintupleVenn(apply(relativePositionResults <=0.3, 2, which),
                       sum(apply(relativePositionResults > 0.3, 1, all, na.rm=T)),
                       'benchmarking_overlap_genes_found_relative_max_0dot3',
                       toolColors)

# Plot differences in whether the gene was found within the first 20% of output genes.
benchmarkQuintupleVenn(apply(relativePositionResults <=0.2, 2, which),
                       sum(apply(relativePositionResults > 0.2, 1, all, na.rm=T)),
                       'benchmarking_overlap_genes_found_relative_max_0dot2',
                       toolColors)

# Plot differences in whether the gene was found within the first 10% of output genes.
benchmarkQuintupleVenn(apply(relativePositionResults <=0.1, 2, which),
                       sum(apply(relativePositionResults > 0.1, 1, all, na.rm=T)),
                       'benchmarking_overlap_genes_found_relative_max_0dot1',
                       toolColors)



###
### Phenotype information comparing the tools (using the phenotypes for which the tools ranked first).
###

# Data transformation for upcoming plots.
dataToPlot <- t(phenotypeCountsWhenToolRankedFirst[which(phenotypeTotals > 5),])
dataYMax <- ceiling(max(apply(dataToPlot, 2, sum))/100)*100+2

# Plot absolute phenotype frequencies grouped by best ranking tool.
initializeGraphicsDevice('benchmarking_first_rank_phenotype_frequencies', width=7, height=8)
par(mar=c(4.5,11.5,0.5,0.5))
barplot(dataToPlot, xaxt='n', yaxt='n', space=0, xlim=c(0,dataYMax), horiz=TRUE)
abline(v=seq(0, dataYMax, 5), col="gray85")
abline(v=seq(0, dataYMax, 20), col="gray70")
barplot(dataToPlot, # excludes input phenotypes with only NA
        col=toolColors, las=1, cex.names=0.5, space=0, xlim=c(0,dataYMax),
        main="", xlab="phenotype input frequency", add=TRUE, horiz=TRUE)
title(ylab="phenotype", mgp=c(10,1,0))
legend("topright", colnames(phenotypeCountsWhenToolRankedFirst), bg="white",
       fill=toolColors)
dev.off()

# Plot relative phenotype frequencies grouped by best ranking tool.
initializeGraphicsDevice('benchmarking_first_rank_phenotype_frequencies_relative', width=7, height=8)
par(mar=c(5.5,11.5,0.0,1.0))
barplot(prop.table(dataToPlot, 2), # excludes input phenotypes with only NA
        col=toolColors, las=1, cex.names=0.5, space=0, main="",
        xlab="phenotype input frequency", border="white", horiz=TRUE)
title(ylab="phenotype", mgp=c(10,1,0))
par(xpd=TRUE) # no clipping for drawing outside plot
legend(-0.38,-1, colnames(phenotypeCountsWhenToolRankedFirst),
       fill=toolColors)
dev.off()

# Ensures dataToPlot cannot be used afterwards.
rm(dataToPlot)



###
### Scatterplots comparing the found positions to the total number of output genes.
###

# Data transformation for upcoming plots.
dataToPlot <- data.frame(hits=unlist(totalResults),
                         position=unlist(positionResults),
                         tool=as.factor(sort(rep(colnames(totalResults), nrow(totalResults)))))

# Plots hits compared to total positions without gado
hitPositionsScatterplot(dataToPlot[which(dataToPlot$tool != "gado"),],
                        'tools_position_totalHits_no_gado',
                        toolColors, 5000, 1000, 2000, 1000,
                        "total number of output genes",
                        "absolute position of validation genes within the output genes",
                        parMar=c(4.5,5.5,0.5,0.5),
                        ylabMgp=c(4,1,0))

# Transforms values to log10 for second plot.
dataToPlot$hits <- log10(dataToPlot$hits)
dataToPlot$position <- log10(dataToPlot$position)

# Plots hits compared to total positions with gado (log10 transformed).
hitPositionsScatterplot(dataToPlot,
                        'tools_position_totalHits_log10', 
                        toolColors, 1, 1, 1, 1,
                        "total number of output genes (log10)",
                        "absolute position of validation genes within the output genes (log10)")

# Ensures dataToPlot cannot be used afterwards.
rm(dataToPlot)



###
### Compares amelie with vibe.
###

# Plots vibe against amelie (absolute).
dataToPlot <- log10(data.frame(vibe=positionResults$vibe,
                               amelie=positionResults$amelie))
plotToolComaprison(dataToPlot,
                   'benchmark_positions_amelie_vibe_absolute',
                   max(dataToPlot, na.rm=TRUE),
                   "absolute position of validation genes in vibe (log10)",
                   "absolute position of validation genes in amelie output genes (log10)")

# Plots vibe against amelie (relative).
dataToPlot <- log10(data.frame(vibe=relativePositionResults$vibe,
                         amelie=relativePositionResults$amelie))
plotToolComaprison(dataToPlot,
                   'benchmark_positions_amelie_vibe_relative',
                   xyMin=floor(min(dataToPlot, na.rm=TRUE)),
                   xyMax=ceiling(max(dataToPlot, na.rm=TRUE)),
                   "relative position of validation genes in vibe (log10)",
                   "relative position of validation genes in amelie (log10)")



###
### Curves showing how many genes were found within different cutoffs.
###

# Plots how often genes were found using cutoffs of total available hits.
# Shows how many hits were found when looking at a specific absolute cutoff
# max positions from all genes found.
plotMatchesFoundWihinRangeCutoff(genePercentageFoundWithinAbsoluteCutoff[1:5000,],
                                 'found_genes_for_absolute_cutoffs',
                                 "absolute cutoff from the output genes",
                                 "fraction of patient cases for which the gene was found",
                                 toolColors, 500, addAbline=FALSE)

# Shows how many hits were found when looking at a specific cutoff fractions
# from all genes found.
plotMatchesFoundWihinRangeCutoff(genePercentageFoundWithinRelativeCutoff,
                                 'found_genes_for_relative_cutoffs',
                                 "relative cutoff from the output genes",
                                 "fraction of patient cases for which the gene was found",
                                 toolColors, 0.1)



###
### Heatmap showing the relative positions from the different tools.
###

# Prepares data for heatmap.
naValue <- 2
dataToPlot <- relativePositionResults
dataToPlot[is.na(dataToPlot)] <- naValue
xClust <- hclust(dist(t(dataToPlot), method = "euclidean"), method="complete")
yClust <- hclust(dist(dataToPlot, method = "euclidean"), method="complete")

# Prepares data for heatmap.
dataToPlot <- cbind(id=rownames(dataToPlot), dataToPlot)
dataToPlot <- melt(dataToPlot, id.vars="id",
                   variable.name="tool",
                   value.name="relativePosition")

# Orders data based on clustering.
dataToPlot$id <- factor(dataToPlot$id,
                        levels = dataToPlot$id[yClust$order],
                        ordered = TRUE)
dataToPlot$tool <- factor(dataToPlot$tool,
                        levels = levels(dataToPlot$tool)[xClust$order],
                        ordered = TRUE)

# Theme for dendrograms.
dendroTheme <- theme(axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.x=element_blank(),
                     axis.text.y=element_blank())

# Tool dendogram. Remove "+ dendroTheme" for comparing values with hmPlot.
xClustPlot <- ggdendrogram(xClust) + dendroTheme

# Patient cases dendogram. Remove "+ dendroTheme" for comparing values with hmPlot.
yClustPlot <- ggdendrogram(yClust, rotate=T) + dendroTheme

# Heatmap.
hmPlot <- ggplot(dataToPlot, aes(x=tool, y=id, colour="")) + # colour is to trick ggplot for "not found" in legend
  theme_bw() +
  theme(legend.box="horizontal",
        legend.box.spacing=unit(0, "cm"),
        legend.key.size=unit(0.4, "cm"),
        legend.text=element_text(size=6),
        legend.title=element_text(size=9),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.text.x=element_text(margin=margin(b=5)),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), # Can be commented to compare with yClustPlot
        axis.ticks.y=element_blank()) +
  labs(y="patient cases",
       fill="relative position",
       colour="not found") +
  geom_tile(aes(fill=relativePosition)) +
  scale_fill_viridis(direction=-1, option="C", na.value="black",
                     limits = c(0,1)) + # excludes the NA value (which is > 1)
  scale_colour_manual(values=NA, labels = "2.00") + # makes sure aes colour is not visible
  scale_y_discrete(expand=c(0,0)) # moves x labels closer

# Filter legend from heatmap and generate plot from it.
hmInfo <- ggplotLegend(hmPlot)
hmInfoPlot <- ggdraw() +
  draw_grob(hmInfo$grobs[[1]], -0.3, 0) +
  draw_grob(hmInfo$grobs[[2]], 0.3, 0)

# Remove legend from original heatmap.
hmPlot <- hmPlot + theme(legend.position = "none")

# Generate image.
initializeGraphicsDevice('heatmap_relative_positions', width=7, height=6)
plot_grid(
  xClustPlot, hmInfoPlot,
  hmPlot, yClustPlot,
  nrow = 2, ncol = 2,
  align="hv",
  axis="tblr",
  rel_widths = c(3,1),
  rel_heights = c(1,2.5),
  scale=c(0.85,1,1,1.08)
)
dev.off()

# Alternative method for heatmap (has export issues + no legend for missing values).
#dataToPlot <- relativePositionResults
#dataToPlot[is.na(dataToPlot)] <- 2
#heatmaply(dataToPlot,
#          col=plasma(100, direction=-1),
#          limits=c(0,1),
#          na.value="black",
#          column_text_angle=0,
#          dist_method="euclidean",
#          hclust_method="complete",
#          xlab="tools",
#          ylab="patient cases",
#          showticklabels=c(TRUE, FALSE),
#          file=paste0(imgExportDir, 'heatmap_relative_positions_heatmaply.png'))



###
### Density plot of z-scores.
###

toolScoresHistogram(toolGeneScores$amelie, "density_amelie", "AMELIE score", c(0,100), 0.04)
toolScoresHistogram(toolGeneScores$gado, "density_gado", "z-score", c(-5,20), 0.25)
toolScoresHistogram(toolGeneScores$vibe, "density_vibe", "gda_max score", c(0,1), 4)
