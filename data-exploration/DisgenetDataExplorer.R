########
# Name:
# DisgenetDataExplorer
#
# Description:
# Creation of plots to explore the DisGeNET dataset.
########

##################
### Libraries  ###
##################

library(plyr)
library(dplyr)
library(RColorBrewer)
library(plotrix)

##################
### Functions  ###
##################

########
# Name:
# sumCounts
#
# Description:
# Sums up the 'countPmid' and 'countNA' columns from a data frame and adds this
# as a new column called 'countTotal'
#
# Input:
# data -  A data frame that has the columns 'countPmid' and 'countNA'
#
# Output:
# A data frame
########
sumCounts <- function(data) {
  data$countTotal <- data$countPmid + data$countNA
  return(data)
}

########
# Name:
# plotPmids
#
# Description:
# Creates a plot that shows how many associations per source have a pmid and
# how many do not (have NA in that field).
#
# Input:
# plotTitle - The title of the plot.
# data      - A data frame with the data to plot. Should contain a 'originalSource'
#             column with the source names and the respective 'countPmid'/'countNA'
#             columns with the actual data that needs to be plotted.
#
# Output:
# N.A. (a plot is generated)
########
plotPmids <- function(plotTitle, data) {
  par(mar=c(7, 4, 4, 2) + 0.1)
  dataToPlot <- log10(t(data[,c('countPmid', 'countNA')]))
  xLabels <- data$source
  barColors <- c('steelblue4', 'slategray2')
  barplot(dataToPlot, beside=T, ylim=c(0,ceiling(max(dataToPlot))), names.arg=xLabels,
          las=2, col=barColors, space=rep(c(1.2,0), 12), main=plotTitle,
          ylab='number of associations (log10)')
  mtext('source', side=1, line=6)
  legend(28,ceiling(max(dataToPlot)), c('with pmid', 'without pmid'), fill=barColors, bty = "n")
}

########
# Name:
# addLevels
#
# Description:
# Adds the levels corresponding to the sources. While for efficiency sourceLevels
# could have been defined outside the function, as it is an integral part of the
# rest of the function it is stored inside it.
#
# Input:
# data  - A data frame containing a column 'source' with sources as defined in
#         the sourceLevels variable.
#
# Output:
# A data frame with an added 'level' column.
########
addLevels <- function(data) {
  # Creates the levels as defined by http://www.disgenet.org/web/DisGeNET/menu/dbinfo#sources
  sourceLevels <- list()
  sourceLevels$curated <- c('UNIPROT', 'CTD_human', 'CLINVAR', 'ORPHANET', 'GWAS', 'PSYGENET',
                            'HPO')
  sourceLevels$animal_models <- c('CTD_mouse', 'CTD_rat', 'MGD', 'RGD')
  sourceLevels$literature <- c('GAD', 'LHGDN', 'BEFREE')
  sourceLevels$order <- c('curated', 'animal models', 'literature')
  
  data$level <- factor(case_when(countsPerSource$source %in% sourceLevels$curated ~ 'curated',
                                 countsPerSource$source %in% sourceLevels$animal_models ~ 'animal models',
                                 countsPerSource$source %in% sourceLevels$literature ~ 'literature'),
                       levels=sourceLevels$order)
  return(data)
}

########
# Name:
# mergeSourcesWithLowCounts
#
# Description:
# Merges rows per level that have a value lower than the cutoff.
#
# Input:
# data    - A data frame with the columns 'source', 'countPmid', 'countNA',
#           'countTotal' and 'level'.
# cutoff  - The cutoff value to determine when rows should be merged (< cutoff).
#
# Output:
# A data frame with merged rows (per level) that had a value lower than cutoff.
########
mergeSourcesWithLowCounts <- function(data, cutoff) {
  # Merges the sources that have less than <cutoff> associations.
  mergedData <- ddply(data[data$countTotal < cutoff,], .(level), summarize,
                      countPmid=sum(countPmid), countNA=sum(countNA),
                      countTotal=sum(countTotal))
  mergedData <- cbind(source=paste('other', mergedData$level), mergedData)
  mergedData <- rbind(data[data$countTotal >= cutoff,], mergedData)
  
  return(mergedData)
}

########
# Name:
# plotAssociationsPerSourceAndLevel
#
# Description:
# Generates a circle diagram with the associations per source and per level.
#
# Input:
# plotTitle         - The title of the plot.
# sourceData        - A data frame with plot data of the sources. Requires a 'source'
#                     and 'countTotal' column.
# levelData         - A data frame with plot data of the levels. Requires a 'level'
#                     and 'countTotal' column.
# outerRadius       - Allows for adjustment of the source (outer) number radii.
# innerRadius       - Allows for adjustment of the source (inner) number radii.
# outerAngles       - Selects which source (outer) values need an angle adjustment.
# outerAngleAdjust  - The angle adjustments for the selected sources.
# innerAngles       - Selects which level (inner) values need an angle adjustment.
# innerAngleAdjust  - The angle adjustments for the selected levels.
#
# Output:
# N.A. (a plot is generated)
########
plotAssociationsPerSourceAndLevel <- function(plotTitle, sourceData, levelData,
                                              outerRadius=0.8, innerRadius=0.3,
                                              outerAngles=c(), outerAngleAdjust=0,
                                              innerAngles=c(), innerAngleAdjust=0) {
  # Generates colors.
  sourceColors <- brewer.pal(nrow(sourceData), 'Set3')
  levelColors <- brewer.pal(nrow(levelData), 'Dark2')
  
  # Create empty pie chart.
  pie(1, col=c('white'), border = NA, labels='', main=plotTitle)
  
  # Adds the 2 layers to the pie chart.
  outerAnglesValues <- floating.pie(-0.6,0, sourceData$countTotal, radius=1,
                                    border=NA, startpos = pi/2, col=sourceColors)
  innerAngleValues <- floating.pie(-0.6,0, levelData$countTotal, radius=0.5,
                                   border=NA, startpos = pi/2, col=levelColors)
  
  # Creates readable slice size numbers.
  outerLabelValues <- format(sourceData$countTotal, big.mark = '.',
                             decimal.mark=',', trim=T)
  innerLabelValues <- format(levelData$countTotal, big.mark = '.',
                             decimal.mark=',', trim=T)
  
  # Optional angle adjustments.
  outerAnglesValues[outerAngles] <- outerAnglesValues[outerAngles] + outerAngleAdjust
  innerAngleValues[innerAngles] <- innerAngleValues[innerAngles] + innerAngleAdjust
  
  # Adds numbers to pie chart.
  pie.labels(-0.6,0, outerAnglesValues, outerLabelValues, radius=outerRadius)
  pie.labels(-0.6,0, innerAngleValues, innerLabelValues, radius=innerRadius)
  
  # Adds legends to pie chart.
  text(0.8,0.3, 'Sources', pos=4, font=2)
  legend(0.8,0.3, sourceData$source[sourceData$countTotal > 0], fill=sourceColors, border=NA, bty = "n")
  text(0.8,0.9, 'Levels', pos=4, font=2)
  legend(0.8,0.9, levelData$level[levelData$countTotal > 0], fill=levelColors, border=NA, bty = "n")
  
  # Calculates, formats and adds total to the pie chart.
  total <- format(sum(sourceData$countTotal), big.mark = '.',
                  decimal.mark=',', trim=T)
  mtext(paste('Total:', total), side=1, at=-0.6)
}


##################
###    Code    ###
##################

# Defaults
oldPar <- par()

######## 
######## Initial data loading/processing.
######## 

# Adjust to local system.
baseDir <- '~/programming/projects/vibe/data/DisGeNET/'
imgExportDir <- '~/Documents/afstuderen/verslag/img/'

# For reading in files:
# Disabled quoting as otherwise not all lines are read ('EOF within quoted string').
# Disabled comment character as some sentences could contain a #, making a row missing columns.

# Loads gene disease pmid association file and creates second variable for phenotype-only associations.
geneDiseasePmidAssociations <- read.table(gzfile(paste0(baseDir, 'all_gene_disease_pmid_associations.tsv.gz')),
                                          header=T, sep='\t', quote="", comment.char="")
geneDiseasePmidPhenotypeAssociations <- geneDiseasePmidAssociations[geneDiseasePmidAssociations$diseaseType == 'phenotype',]

# Merges the counts per source in non-NA's and NA's + adds a total column.
countsPerSource <- ddply(geneDiseasePmidAssociations, .(source=originalSource), summarize,
                         countPmid=sum(!is.na(pmid)), countNA=sum(is.na(pmid)))
countsPerSource <- sumCounts(countsPerSource)
countsPerSource.phenotype <- ddply(geneDiseasePmidPhenotypeAssociations, .(source=originalSource), summarize,
                                   countPmid=sum(!is.na(pmid)), countNA=sum(is.na(pmid)))
countsPerSource.phenotype <- sumCounts(countsPerSource.phenotype)

# Adds the levels to the data.
countsPerSource <- addLevels(countsPerSource)
countsPerSource.phenotype <- addLevels(countsPerSource.phenotype)




######## 
######## Plots the number of non-NA and NA associations per source.
########
pdf(paste0(imgExportDir, 'pubmed-count-per-source-all.pdf'), 7,5)
plotPmids('Total number of gene associations per source', countsPerSource)
dev.off()

pdf(paste0(imgExportDir, 'pubmed-count-per-source-phenotype.pdf'), 7,5)
plotPmids('Total number of gene-phenotype associations per source', countsPerSource.phenotype)
dev.off()

######## 
######## Plotting the total number of phenotype gene-disease associations per source/level.
######## 
# Calculate number of associations per level.
countsPerLevel <- ddply(countsPerSource, "level", numcolwise(sum))
countsPerLevel.phenotype <- ddply(countsPerSource.phenotype, "level", numcolwise(sum))

# Merges the sources that have less than 1000 associations.
MergedCountsPerSource <- mergeSourcesWithLowCounts(countsPerSource, 20000)
MergedCountsPerSource.phenotype <- mergeSourcesWithLowCounts(countsPerSource.phenotype, 200)

# Sort the sources on level.
MergedCountsPerSource <- arrange(MergedCountsPerSource, level)
MergedCountsPerSource.phenotype <- arrange(MergedCountsPerSource.phenotype, level)

# Generates plots.
pdf(paste0(imgExportDir, 'gene-associations.pdf'), 9,7)
plotAssociationsPerSourceAndLevel('number of gene associations within DisGeNET 5.0',
                                  MergedCountsPerSource, countsPerLevel,
                                  outerRadius=c(0.9, 0.84, 0.76, 0.68, 0.8, 0.7, 0.79),
                                  innerRadius=c(0.35,0.25,0.15),
                                  outerAngles=c(3,4,6),
                                  outerAngleAdjust=c(-0.07,0.02,0.05),
                                  innerAngles=1,
                                  innerAngleAdjust=-0.08)
dev.off()

pdf(paste0(imgExportDir, 'gene-phenotype-associations.pdf'), 9,7)
plotAssociationsPerSourceAndLevel('number of gene-phenotype associations within DisGeNET 5.0',
                                  MergedCountsPerSource.phenotype, countsPerLevel.phenotype,
                                  outerRadius=c(0.9, 0.65, 0.85, 0.7, 0.55, 0.7, 0.7, 0.75),
                                  innerRadius=c(0.2,0.3,0.15),
                                  outerAngles=c(3,4,5),
                                  outerAngleAdjust=c(-0.05,0.02,0.1))
dev.off()
