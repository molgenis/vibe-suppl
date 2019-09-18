########
# Name:
# ToolComparisonResultsGenerator.R
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
imgExportDir <- paste0(baseImgExportDir, "tools/")

# Defines color palette.
toolColors <- brewer.pal(5, 'Set2') # Adjust if changing number of tools!!!

# Load data.
amelie <- readResultFile("results/amelie.tsv")
gado <- readResultFile("results/gene_network.tsv")
phenomizer <- readResultFile("results/phenomizer.tsv")
phenotips <- readResultFile("results/phenotips.tsv")
vibe <- readResultFile("results/vibe_2018-07-06_gda_max.tsv")

# Sorts benchmark results so that row order is identical.
amelie <- sortRows(amelie)
gado <- sortRows(gado)
phenomizer <- sortRows(phenomizer)
phenotips <- sortRows(phenotips)
vibe <- sortRows(vibe)

# Loads in best score per tool.
toolGeneScores <- read.table(paste0(baseDir, "results/tool_gene_scores.tsv"),
                             header=T, sep="\t")
toolGeneScores <- toolGeneScores[,3:5]



###
### Calculations for comparison of different tools (with best scoring vibe result) - positions
###

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



###
### Calculations for comparison of different tools (with best scoring vibe result) - rankings
###

# Calculate absolute tool rankings.
toolRankings <- calculateRankings(positionResults)
toolRankingCounts <- calculateRankingCounts(positionResults, toolRankings)



###
### Calculations for comparison of different tools (with best scoring vibe result) - phenotypes
###

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
rm(foundCounts, toolName)

# The frequencies of the input phenotypes for which a tool ranked first (no all NA).
phenotypeTotals <- apply(phenotypeCountsWhenToolRankedFirst, 1, sum)



###
### Calculations for comparison of different tools (with best scoring vibe result) - found matches within position cutoffs
###

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
rm(x, highestPossibleResultForTool)

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

rm(naExceptPhenomizerOnly)



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
rm(dataToPlot, dataYMax)



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

rm(dataToPlot)



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

rm(naValue, dataToPlot, xClust, yClust, xClustPlot, yClustPlot, hmPlot, hmInfo, hmInfoPlot, dendroTheme)

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
toolScoresHistogram(toolGeneScores$gado, "density_gado", "z-score", c(-5.5,20), 0.25)
toolScoresHistogram(toolGeneScores$vibe, "density_vibe", "gda_max score", c(0,1), 4)

###
### Removes any script-sepcific variables.
###

rm(toolColors, imgExportDir, amelie, gado, phenomizer, phenotips, vibe,
   genePercentageFoundWithinAbsoluteCutoff, genePercentageFoundWithinRelativeCutoff,
   genesFoundWithinAbsoluteCutoff, genesFoundWithinRelativeCutoff,
   allPhenotypeNames, phenotypeCountsWhenToolRankedFirst, phenotypeTotals,
   positionResults, relativePositionResults,
   toolGeneScores,
   toolRankings, toolRankingCounts,
   totalResults)
