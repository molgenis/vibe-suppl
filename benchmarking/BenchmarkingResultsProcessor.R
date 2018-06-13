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
        strsplit(benchmarkResults[benchmarkDataRow["lovd"],], ",")[[1]])
}



##################
###    Code    ###
##################

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

phenotips <- readResultFile("results/phenotips.tsv")
phenomizer <- readResultFile("results/phenomizer.tsv")
amelie <- readResultFile("results/amelie.tsv")
vibe.20180503 <- readResultFile("results/vibe_2018-05-02.tsv")

# Process data
positionResults <- data.frame(amelie=resultsPositionCalculator(benchmarkData, amelie),
                              phenomizer=resultsPositionCalculator(benchmarkData, phenomizer),
                              phenotips=resultsPositionCalculator(benchmarkData, phenotips),
                              vibe.20180503=resultsPositionCalculator(benchmarkData, vibe.20180503))

# Count missing values.
resultsWithNa <- positionResults[apply(apply(positionResults, c(1,2), is.na), 1, any),]
naCounts <- apply(apply(resultsWithNa, c(1,2), is.na), 2, sum)

### Relative position calculations.
positionOrderedToolNames <- sapply(apply(positionResults, 1, sort, na.last=NA), names)
relativePositionResults <- sapply(1:4, function(x) { table(sapply(positionOrderedToolNames, '[', x)) })
colnames(relativePositionResults) <- 1:4
relativePositionResults <- cbind(relativePositionResults, naCounts)
apply(relativePositionResults, 1, sum)

###
### Boxplot comparing absolute positions of genes.
###
yAxisMax <- ceiling(max(positionResults, na.rm=T)/100)*100

postscript(paste0(imgExportDir, 'benchmarking_gene_position_boxplots.eps'), width=6, height=8)
boxplot(positionResults, xaxt='n',yaxt='n', pch="")
abline(h=1, col="gray92")
abline(h=seq(100, yAxisMax, 100), col="gray92")
boxplot(positionResults, las=1, pch=20, yaxt='n',
        main="position of relevant genes among different tools",
        xlab="tool", ylab="position", col="white", add=T)
axis(2, las=1, at=1)
axis(2, las=1, at=seq(500, yAxisMax,500))
dev.off()

###
### Barplot for relative positions.
###

barplot(t(relativePositionResults), col=rev(brewer.pal(5, 'RdYlGn')))

###
### Venn diagrams showing differences between tools.
###
colors <- brewer.pal(ncol(positionResults), 'Set3')
#grid.newpage()

# Plot differences in whether the gene was found.
cairo_ps(paste0(imgExportDir, 'benchmarking_genes_found.eps'), width=8, height=8)
grid.draw(venn.diagram(apply(!is.na(positionResults), 2, which), NULL,
                       main.fontfamily="Helvetica", sub.fontfamily="Helvetica",
                       fontfamily="Helvetica", cat.fontfamily="Helvetica",
                       fill=colors))
dev.off()

# Plot differences in whether the gene was found within the first 100 positions.
cairo_ps(paste0(imgExportDir, 'benchmarking_genes_max_100.eps'), width=8, height=8)
grid.draw(venn.diagram(apply(positionResults <=100, 2, which), NULL,
                       main.fontfamily="Helvetica", sub.fontfamily="Helvetica",
                       fontfamily="Helvetica", cat.fontfamily="Helvetica",
                       fill=colors))
dev.off()

# Plot differences in whether the gene was found within the first 20 positions.
cairo_ps(paste0(imgExportDir, 'benchmarking_genes_max_020.eps'), width=8, height=8)
grid.draw(venn.diagram(apply(positionResults <=20, 2, which), NULL,
                       main.fontfamily="Helvetica", sub.fontfamily="Helvetica",
                       fontfamily="Helvetica", cat.fontfamily="Helvetica",
                       fill=colors))
dev.off()
