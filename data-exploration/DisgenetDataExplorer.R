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
# calculateCountsPerSource
#
# Description:
# Calculates the count per source.
#
# Input:
# data - a data frame with for each association a row, a column
#        'originalSource' describing the source for each row and a column 'pmid'
#        containing either the PMID of the source or NA if not available.
#
# Output:
# A data frame with the different sources as rows and the tool names, PMID counts,
# NA counts (when PMID was NA) and a total of these two as columns.
########
calculateCountsPerSource <- function(data) {
  countsPerSource <- ddply(data, .(source=originalSource), summarize,
                           countPmid=sum(!is.na(pmid)), countNA=sum(is.na(pmid)))
  countsPerSource <- sumCounts(countsPerSource)
  countsPerSource <- addLevels(countsPerSource)
  return(countsPerSource)
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
  par(mar=c(7, 4, 1, 0) + 0.1)
  dataToPlot <- log10(t(data[,c('countPmid', 'countNA')]))
  xLabels <- data$source
  barColors <- c('steelblue4', 'slategray2')
  barplot(dataToPlot, beside=T, ylim=c(0,ceiling(max(dataToPlot))), names.arg=xLabels,
          las=2, col=barColors, space=rep(c(1.2,0), 12), main=plotTitle,
          ylab='number of associations (log10)')
  mtext('source', side=1, line=6)
  legend(28,ceiling(max(dataToPlot)), c('with pmid', 'without pmid'), fill=barColors, bty = "n")
  par(oldPar)
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
  
  data$level <- factor(case_when(data$source %in% sourceLevels$curated ~ 'curated',
                                 data$source %in% sourceLevels$animal_models ~ 'animal models',
                                 data$source %in% sourceLevels$literature ~ 'literature'),
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
  
  par(mar=c(1.5, 0, 0, 0) + 0.1)
  # Create empty pie chart.
  pie(1, col=c('white'), border = NA, labels='', main=plotTitle)
  
  # Adds the 2 layers to the pie chart.
  outerAngleValues <- floating.pie(-0.6,0, sourceData$countTotal, radius=1,
                                    border=NA, startpos = pi/2, col=sourceColors)
  innerAngleValues <- floating.pie(-0.6,0, levelData$countTotal, radius=0.5,
                                   border=NA, startpos = pi/2, col=levelColors)
  
  # Creates readable slice size numbers.
  outerLabelValues <- format(sourceData$countTotal, big.mark = '.',
                             decimal.mark=',', trim=T)
  innerLabelValues <- format(levelData$countTotal, big.mark = '.',
                             decimal.mark=',', trim=T)
  
  # Optional angle adjustments.
  outerAngleValues[outerAngles] <- outerAngleValues[outerAngles] + outerAngleAdjust
  innerAngleValues[innerAngles] <- innerAngleValues[innerAngles] + innerAngleAdjust
  
  # Adds numbers to pie chart.
  pie.labels(-0.6,0, outerAngleValues, outerLabelValues, radius=outerRadius)
  pie.labels(-0.6,0, innerAngleValues, innerLabelValues, radius=innerRadius)
  
  # Adds legends to pie chart.
  text(0.8,0.3, 'Sources', pos=4, font=2)
  legend(0.8,0.3, sourceData$source[sourceData$countTotal > 0], fill=sourceColors, bty = "n")
  text(0.8,0.9, 'Levels', pos=4, font=2)
  legend(0.8,0.9, levelData$level[levelData$countTotal > 0], fill=levelColors, bty = "n")
  
  # Calculates, formats and adds total to the pie chart.
  total <- format(sum(sourceData$countTotal), big.mark = '.',
                  decimal.mark=',', trim=T)
  mtext(paste('Total:', total), side=1, at=-0.6)
  par(oldPar)
}

########
# Name:
# mergeValues
#
# Description:
# Sums values in the defined slice sizes. If the last slice is smaller than the
# slice size, these values are still treated as a single slice.
#
# Input:
# data              - A vector containing integers.
# sliceSize         - A single integer defining the size of the slices.
#
# Output:
# A vector containing integer with labels describing the range each value is
# created from.
########
mergeValues <- function(data, sliceSize) {
  dataLength <- length(data)
  fullSliceLimit <- dataLength - dataLength %% sliceSize
  
  resultData <- apply(matrix(data[1:fullSliceLimit], byrow=T, ncol=sliceSize), 1, sum)
  resultLabels <- paste(seq(1, fullSliceLimit, sliceSize),
                           seq(sliceSize, fullSliceLimit, sliceSize), sep='-')
  
  if(fullSliceLimit < dataLength) {
    restData <- sum(data[(fullSliceLimit+1):dataLength])
    resultData <- c(resultData,restData)
    
    restLabel <- paste(fullSliceLimit+1, fullSliceLimit+sliceSize, sep='-')
    names(resultData) <- c(resultLabels, restLabel)
  }
  
  return(resultData)
}

##################
###    Code    ###
##################

# Defaults
oldPar <- par()
setEPS() # Sets EPS engine for writing images.

######## 
######## Initial data loading/processing.
######## 

# Adjust to local system.
baseDir <- '~/programming/projects/vibe/data/DisGeNET/'
imgExportDir <- '~/Documents/afstuderen/verslag/img/'

# For reading in files:
# Disabled quoting as otherwise not all lines are read ('EOF within quoted string').
# Disabled comment character as some sentences could contain a #, making a row missing columns.

# Loads gene disease pmid association file and creates subsets.
geneAssociationsPmid <- read.table(gzfile(paste0(baseDir, 'all_gene_disease_pmid_associations.tsv.gz')),
                               header=T, sep='\t', quote="", comment.char="")
geneDiseaseAssociationsPmid <- geneAssociationsPmid[geneAssociationsPmid$diseaseType == 'disease',]
geneGroupAssociationsPmid <- geneAssociationsPmid[geneAssociationsPmid$diseaseType == 'group',]
genePhenotypeAssociationsPmid <- geneAssociationsPmid[geneAssociationsPmid$diseaseType == 'phenotype',]

# Calculates the counts per source in non-NA's and NA's PMIDs + adds a total column.
geneAssociationsSourceCounts <- calculateCountsPerSource(geneAssociationsPmid)
geneDiseaseAssociationsSourceCounts <- calculateCountsPerSource(geneDiseaseAssociationsPmid)
geneGroupAssociationsSourceCounts <- calculateCountsPerSource(geneGroupAssociationsPmid)
genePhenotypeAssociationsSourceCounts <- calculateCountsPerSource(genePhenotypeAssociationsPmid)

# Loads variant disease pmid association file and creates second variable for phenotype-only associations.
variantAssociationsPmid <- read.table(gzfile(paste0(baseDir, 'all_variant_disease_pmid_associations.tsv.gz')),
                                             header=T, sep='\t', quote="", comment.char="")
variantDiseaseAssociationsPmid <- variantAssociationsPmid[variantAssociationsPmid$diseaseType == 'disease',]
variantGroupAssociationsPmid <- variantAssociationsPmid[variantAssociationsPmid$diseaseType == 'group',]
variantPhenotypeAssociationsPmid <- variantAssociationsPmid[variantAssociationsPmid$diseaseType == 'phenotype',]

# Loads association files created using bash from the "all_gene_disease_pmid_associations.tsv.gz" file.
geneAssociations <- read.table(gzfile(paste0(baseDir, 'complete_gene_associations.tsv.gz')),
                               header=T, sep='\t', quote="", comment.char="")



######## 
######## Shows the number of (phenotype) gene/variant pmid associations.
########
nrow(geneDiseaseAssociationsPmid)
nrow(geneGroupAssociationsPmid)
nrow(genePhenotypeAssociationsPmid)
nrow(geneAssociationsPmid)

nrow(variantDiseaseAssociationsPmid)
nrow(variantGroupAssociationsPmid)
nrow(variantPhenotypeAssociationsPmid)
nrow(variantAssociationsPmid)

######## 
######## Shows the possible disease types.
########
unique(geneAssociationsPmid$diseaseType)
unique(variantAssociationsPmid$diseaseType)

######## 
######## Shows number of unique genes for the associations.
########
dataToPlot <- c(all=length(unique(geneAssociationsPmid$geneId)),
                disease=length(unique(geneDiseaseAssociationsPmid$geneId)),
                group=length(unique(geneGroupAssociationsPmid$geneId)),
                phenotype=length(unique(genePhenotypeAssociationsPmid$geneId)))

ylimTop <- ceiling(max(dataToPlot/1000))
postscript(paste0(imgExportDir, 'unique_genes.eps'), width=5, height=5)
barplot(dataToPlot/1000, ylim=c(0,20), las=1,
        ylab='number of unique genes (per thousand)', col="steelblue4")
dev.off()

rm(dataToPlot, ylimTop)

######## 
######## Plots the number of non-NA and NA associations per source.
########
postscript(paste0(imgExportDir, 'pubmed_count_per_source_all.eps'), width=7, height=4)
plotPmids('', geneAssociationsSourceCounts)
dev.off()

######## 
######## Plotting the total number of phenotype gene-disease associations per source/level.
######## 
postscript(paste0(imgExportDir, 'gene_associations.eps'), width=8, height=4.5)
plotAssociationsPerSourceAndLevel('',
                                  # Counts per source (sorted and sources with <2000 are merged).
                                  arrange(mergeSourcesWithLowCounts(geneAssociationsSourceCounts, 20000), level),
                                  # Counts per level.
                                  ddply(geneAssociationsSourceCounts, "level", numcolwise(sum)),
                                  outerRadius=c(0.9, 0.84, 0.76, 0.68, 0.8, 0.7, 0.79),
                                  innerRadius=c(0.35,0.25,0.15),
                                  outerAngles=c(3,4,6),
                                  outerAngleAdjust=c(-0.07,0.02,0.05),
                                  innerAngles=1,
                                  innerAngleAdjust=-0.08)
dev.off()

postscript(paste0(imgExportDir, 'gene_disease_associations.eps'), width=8, height=4.5)
plotAssociationsPerSourceAndLevel('',
                                  # Counts per source (sorted and sources with <2000 are merged).
                                  arrange(mergeSourcesWithLowCounts(geneDiseaseAssociationsSourceCounts, 20000), level),
                                  # Counts per level.
                                  ddply(geneDiseaseAssociationsSourceCounts, "level", numcolwise(sum)),
                                  outerRadius=c(0.9, 0.8, 0.7, 0.6, 0.8, 0.8, 0.6),
                                  innerRadius=c(0.35,0.25,0.15),
                                  outerAngles=c(3,4),
                                  outerAngleAdjust=c(-0.07,0.02))
dev.off()

postscript(paste0(imgExportDir, 'gene_group_associations.eps'), width=8, height=4.5)
plotAssociationsPerSourceAndLevel('',
                                  # Counts per source (sorted and sources with <2000 are merged).
                                  arrange(mergeSourcesWithLowCounts(geneGroupAssociationsSourceCounts, 5000), level),
                                  # Counts per level.
                                  ddply(geneGroupAssociationsSourceCounts, "level", numcolwise(sum)),
                                  outerRadius=c(0.86, 0.78, 0.68, 0.7, 0.7, 0.8),
                                  innerRadius=c(0.35,0.25,0.15))
dev.off()

postscript(paste0(imgExportDir, 'gene_phenotype_associations.eps'), width=8, height=4.5)
plotAssociationsPerSourceAndLevel('',
                                  # Counts per source (sorted and sources with <2000 are merged).
                                  arrange(mergeSourcesWithLowCounts(genePhenotypeAssociationsSourceCounts, 200), level),
                                  # Counts per level.
                                  ddply(genePhenotypeAssociationsSourceCounts, "level", numcolwise(sum)),
                                  outerRadius=c(0.9, 0.65, 0.85, 0.68, 0.55, 0.7, 0.75, 0.75),
                                  innerRadius=c(0.2,0.3,0.15),
                                  outerAngles=c(3,4,5),
                                  outerAngleAdjust=c(-0.05,0.02,0.1))
dev.off()

######## 
######## Plots the number of associations per unique gene-disease combination.
########
associationCounts <- tabulate(geneAssociations$nAssociations)

# If no 'BEFREE=' present, resulting output cannot be cast to an integer.
# This means each NA is gene-disease combination WITHOUT a befree association.
befreeAssociationCounts <- as.integer(gsub(";.*", "", gsub(".*BEFREE=", "", geneAssociations$sources)))
befreeAssociationCounts[which(is.na(befreeAssociationCounts))] <- 0 # Sets NAs to 0
associationCountsWithoutBefree <- tabulate(geneAssociations$nAssociations - befreeAssociationCounts)

barColors <- c('steelblue4', 'slategray2')

postscript(paste0(imgExportDir, 'associations_per_unique_gene_disease_combination.eps'), width=8, height=4.5)
xLim <- 200
barplot(log10(associationCounts[1:xLim]), ylim=c(0,6),
        border=NA, col=barColors[2], space=0, las=1,
        xlab='number of assocations',
        ylab='number of unique gene-disease combinations (log10)')
barplot(log10(associationCountsWithoutBefree[1:xLim]), ylim=c(0,6), border=NA,
        col=barColors[1], space=0, las=1, add=T, yaxt='n')
xLabels <- c(1,seq(20,xLim,20))
axis(1, at=xLabels-0.5, xLabels)
legend(140,6, c('without BeFree', 'with BeFree'), fill=barColors, bty = "n")
dev.off()

postscript(paste0(imgExportDir, 'associations_per_unique_gene_disease_combination_merged.eps'), width=8, height=4.5)
dataToPlotWithBefree <- log10(mergeValues(associationCounts,200))
dataToPlotWithoutBefree <- log10(mergeValues(associationCountsWithoutBefree,200))
x <- barplot(dataToPlotWithBefree, ylim=c(0,6),
        border=NA, col=barColors[2], las=1, xaxt='n',
        ylab='number of unique gene-disease combinations (log10)')
barplot(dataToPlotWithoutBefree, ylim=c(0,6), border=NA,
        col=barColors[1], las=1, add=T, yaxt='n', xaxt='n')
text(cex=1, x=x+0.5, y=-0.25, names(dataToPlotWithBefree), xpd=TRUE, srt=45, pos=2)
mtext(side=1, line=4, 'number of assocations')
legend(14,6, c('without BeFree', 'with BeFree'), fill=barColors, bty = "n")
dev.off()

rm(x, associationCounts, befreeAssociationCounts, associationCountsWithoutBefree, barColors, xLim)
