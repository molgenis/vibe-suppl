########
# Name:
# PaperPlots.R
#
# Description:
# Generates plots for the paper from the result data. This script functions
# completely independent and does not require the "BenchmarkResultsGenerics.R"
# script.
# 
# Important:
# Be sure to adjust the items in the config section before running the script!
########

##################
### Config     ###
##################

# Directory storing script input:
# - benchmark_data-ncbi_id.tsv
# - amelie.tsv
# - gado.tsv
# - hiphive.tsv
# - phenix.tsv
# - phenomizer.tsv
# - phenotips.tsv
# - pubcasefinder.tsv
# - vibe.tsv
baseDir <- '~/Desktop/zenodo_download/'

# Directory to write the output images to (be sure that directoy exists).
imgExportDir <- '~/Desktop/zenodo_download/out/'





##################
### Libraries  ###
##################

library(dplyr)
library(reshape2)
library(ggplot2)
library(ggdendro)
library(cowplot)





##################
### Functions  ###
##################

########
# Name:
# initializeGraphicsDevice
#
# Description:
# A simple wrapper for writing a file to a predetermined graphics device
# that prepends the directory and appends the file extension.
#
# Input:
# fileName - Filename (excluding file extension) to be used for storage.
# width	- The width of the device in inches. -> see ?pdf()
# height - The height of the device in inches. -> see ?pdf()
#
# Output:
#
########
initializeGraphicsDevice <- function(fileName, width, height) {
  pdf(paste0(imgExportDir, fileName, ".pdf"), width=width, height=height)
}

########
# Name:
# ggSaveCustom
#
# Description:
# A simple wrapper for writing a file using ggsave that prepends the directory
# and appends the file extension.
#
# Input:
# fileName - Filename (excluding file extension) to be used for storage.
# width	- The width of the device in inches. -> see ?ggsave()
# height - The height of the device in inches. -> see ?ggsave()
#
# Output:
#
########
ggSaveCustom <- function(fileName, width, height) {
  ggsave(paste0(imgExportDir, fileName, ".pdf"), width=width, height=height)
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





##################
###    Code    ###
##################

###
### Basic data preperations for further usage.
###

# Load benchmark input data (with ncbi gene id's).
benchmarkData <- read.table(paste0(baseDir,"benchmark_data-ncbi_id.tsv"), header=T,
                            sep="\t",colClasses=c(rep("character", 3),
                                                  "factor", "character"))

# Loads benchmark results and sorts it so that row order is identical.
amelie <- sortRows(readResultFile("amelie.tsv"))
gado <- sortRows(readResultFile("gado.tsv"))
phenomizer <- sortRows(readResultFile("phenomizer.tsv"))
phenotips <- sortRows(readResultFile("phenotips.tsv"))
hiphive <- sortRows(readResultFile("hiphive.tsv"))
phenix <- sortRows(readResultFile("phenix.tsv"))
pubcf <- sortRows(readResultFile("pubcasefinder.tsv"))
vibe <- sortRows(readResultFile("vibe.tsv"))

# Calculate absolute positions. This is done by processing benchmarkData row-by-row
# and therefore the LOVD row order of positionResults is equal to that of benchmarkData.
positionResults <- data.frame(amelie=resultsPositionCalculator(benchmarkData, amelie),
                              gado=resultsPositionCalculator(benchmarkData, gado),
                              phenomizer=resultsPositionCalculator(benchmarkData, phenomizer),
                              phenotips=resultsPositionCalculator(benchmarkData, phenotips),
                              vibe=resultsPositionCalculator(benchmarkData, vibe),
                              hiphive=resultsPositionCalculator(benchmarkData, hiphive),
                              pubcf=resultsPositionCalculator(benchmarkData, pubcf),
                              phenix=resultsPositionCalculator(benchmarkData, phenix))

# Calculate total number of genes.
totalResults <- data.frame(amelie=calculateTotalGenesFound(amelie),
                           gado=calculateTotalGenesFound(gado),
                           phenomizer=calculateTotalGenesFound(phenomizer),
                           phenotips=calculateTotalGenesFound(phenotips),
                           vibe=calculateTotalGenesFound(vibe),
                           hiphive=calculateTotalGenesFound(hiphive),
                           pubcf=calculateTotalGenesFound(pubcf),
                           phenix=calculateTotalGenesFound(phenix),
                           row.names=rownames(amelie))

# Replicates some of the totalResults so that size is equal to positionResults.
totalResults <- totalResults[benchmarkData$lovd,]

###
### Data info.
###

# Check whether the rows of all benchmarks are identically ordered.
rownames(amelie) == rownames(gado) &&
  rownames(gado) == rownames(phenomizer) &&
  rownames(phenomizer) == rownames(phenotips) &&
  rownames(phenotips) == rownames(vibe) &&
  rownames(vibe) == rownames(hiphive) &&
  rownames(hiphive) == rownames(phenix) &&
  rownames(phenix) == rownames(pubcf)

# Check whether there is any phenotype set for which no single gene was found.
any(is.na(totalResults))

###
### Figure 1: Scatterplot with means and missing
###

# Preperations.
colnames(positionResults) <- c("AMELIE", "GADO", "Phenomizer", "Phenotips", "VIBE", "hiPHIVE",  "PubCF.", "PhenIX")
posRelM <- melt(positionResults, id.vars = 0)
totResM <- melt(totalResults, id.vars = 0)
posRelM$total <- totResM$value
colnames(posRelM) <- c("tool", "rank", "total")
posRelM$relative <- posRelM$rank / posRelM$total
levels(posRelM$tool)[match("PubCF.",levels(posRelM$tool))] <- "PubCaseF."

gd <- posRelM %>% 
  group_by(tool) %>% 
  summarise(total = mean(total, na.rm=T),
            rank  = mean(rank, na.rm=T))
toolNaRanks <- aggregate(rank ~ tool, data=posRelM, function(x) {sum(is.na(x))}, na.action = NULL)
gd$NAs <- paste(toolNaRanks$tool, " (", toolNaRanks$rank, " missed)", sep="")

colours <- c("Phenomizer" = "#CC79A7",
             "Phenotips" = "#D55E00",
             "PhenIX" = "#009E73",
             "AMELIE" = "#F0E442",
             "VIBE" = "#0072B2",
             "PubCaseF." = "#56B4E9",
             "hiPHIVE" = "#E69F00",
             "GADO" = "#505050") # http://jfly.uni-koeln.de/color/#pallet
gd$labX <- c(250,   3100,    10,    10,   250,  180,   10,   10)
gd$labY <- c(30000, 30000, 30000, 3000, 10000, 3000, 10000, 1000)

# Plotting figure.
ggplot() +
  geom_point(data = posRelM, aes(x=total, y=rank, color=tool), size = 0.3) + #0.1
  geom_point(data = gd, aes(x=total, y=rank, fill=tool), shape = 21, color = "black", stroke = 1, size = 2) +
  # legend stuff
  geom_point(aes(x=11, y=150), color = "black", size=0.3) +
  geom_text(aes(x=11, y=150, label = "= one causal gene"), color = "black", size = 2, hjust = 0, nudge_x = 0.1) +
  geom_point(aes(x=11, y=300), shape = 1, color = "black", stroke = 1, size = 2) +
  geom_text(aes(x=11, y=300, label="= tool X and Y means"), color="black", hjust = 0, size = 2, nudge_x = 0.1) +
  geom_text(aes(x=8.5, y=30, label="Total: 308\ncausal genes"), color="black", hjust = 0, size = 2, nudge_x = 0.1) +
  geom_label(data = gd, aes(x = labX, y = labY, label = NAs, fill = tool), color="white", hjust = 0, size = 2, fontface = "bold") +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) + scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 40000)) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black"), axis.ticks = element_line(colour = "black"), legend.position = "none", axis.text = element_text(color = "black")) +
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours) +
  labs(x = "Total number of candidate genes returned", y = "Rank of the causal gene")
ggSaveCustom("Figure1", width=4, height=2.5)

###
### Figure 2: Heatmap showing absolute positions from the different tools.
###

# Preperations.
customBreaks <- c(1, 10, 100, 1000, 10000, 40000)
customColours <- c("#0072B2","#009E73","#F0E442", "#D55E00")
naValue <- 50000
dataToPlot <- positionResults
dataToPlot[is.na(dataToPlot)] <- naValue

xClust <- hclust(dist(t(dataToPlot), method = "euclidean"), method="complete")
yClust <- hclust(dist(dataToPlot, method = "euclidean"), method="complete")

dataToPlot <- cbind(id=rownames(dataToPlot), dataToPlot)
dataToPlot <- melt(dataToPlot, id.vars="id",
                   variable.name="tool",
                   value.name="absPosition")

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
                     axis.text.y=element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))

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
        axis.title = element_text(size = 20),
        axis.text.x=element_text(colour="black", margin=margin(b=5)),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), # Can be commented to compare with yClustPlot
        axis.ticks.y=element_blank()) +
  labs(x="Gene prioritization tool", y="Causal genes in patient cases",
       fill="Causal gene rank\n1=best, 40k=worst",
       colour="Missed gene\n(not found)") +
  geom_tile(aes(fill=absPosition)) +
  scale_fill_gradientn(colours = customColours, trans = "log", limit=c(1,40000), breaks = customBreaks, na.value = "black") +
  scale_colour_manual(values=NA, labels = "") + # makes sure aes colour is not visible
  scale_y_discrete(expand=c(0,0)) # moves x labels closer

# Filter legend from heatmap and generate plot from it.
hmInfo <- ggplotLegend(hmPlot)
hmInfoPlot <- ggdraw() +
  draw_grob(hmInfo$grobs[[2]], 0.25, 0) + # "not found" 0.1,0
  draw_grob(hmInfo$grobs[[1]], -0.70, 0) + # color gradient legend
  draw_label("X/Y dendrograms by\nhierarchical clustering\non Euclidean distance", x = 0.51, y = 0.15, size = 8)

# Remove legend from original heatmap.
hmPlot <- hmPlot + theme(legend.position = "none")

# Generate image.
initializeGraphicsDevice('Figure2', width=7, height=6)
plot_grid(
  xClustPlot, hmInfoPlot,
  hmPlot, yClustPlot,
  nrow = 2, ncol = 2,
  align="hv",
  axis="tblr",
  rel_widths = c(3,1),
  rel_heights = c(1,2.5),
  scale=c(0.9,0.95,1,1.08)
)
dev.off()

###
### Unique found per tool at defined cutoff.
###

# Finds all cases found at cutoff. 
cutoffAbs <- 20
withinTopX <- apply(positionResults <= cutoffAbs, 2, which)

# Finds all unique cases at cutoff.
uniqueFound <- sapply(1:length(withinTopX), function(x) {
  withinTopX[[x]][which(!unlist(withinTopX[x]) %in% unlist(withinTopX[-x]))]
})
names(uniqueFound) <- names(withinTopX)

# Number of unique found per tool.
sapply(uniqueFound, length)

# Total unique.
length(unlist(uniqueFound))

# Total solved at cutoff.
length(unique(unlist(withinTopX)))

#############
# FIGURE 3
# extra analysis to show practical value
#############

# load NCBI symbols of CGD genes and convert to vector
cgdDF <- read.table("/Users/joeri/github/vibe-suppl/cgd/CGD_4feb2020-ncbi-gene-ids-minimized.tsv", header = T)
cgd <- cgdDF[,1]

# set seed for pseudo-random numbers for reproducibility
# don't forget to reset seed when re-running code
set.seed(0) 

# empty dataframe to capture results
out <- data.frame()

# the names of the dataframes used as inputs
dfs <- c("gado", "vibe", "hiphive", "phenix", "pubcf", "phenomizer", "phenotips", "amelie")

# iterate over the tools
for(df in dfs)
{
  data <- get(df)
  # iterate over the patient cases in benchmark
  for (i in 1:nrow(benchmarkData)) {
    row <- benchmarkData[i,]
    permGI <- c()
    # do permutations to stabilize results, use odd number to prevent ties in median
    for (p in 1:5) { 
      # rank of the causal gene
      geneRank <- match(row$gene, strsplit(data[unlist(row$lovd),"suggested_genes"], ",")[[1]])
      # if not in tool output, set NA and continue
      if(is.na(geneRank)) {
        geneIndex <- NA
      } else {
        # select all CGD genes minus the current causal one
        cgdMinusCurrentGene <- cgd[-which(cgd==row$gene)]
        # randomly sample 19 other genes
        samp <- sample(cgdMinusCurrentGene, 19)
        # get the ranks of the random genes
        sampRanks <- match(samp, strsplit(data[unlist(row$lovd),"suggested_genes"], ",")[[1]])
        # collate, sort, and get the index of the causal gene, 1 is best, 20 worst
        allRanks <- c(sampRanks, geneRank)
        sortRanks <- sort(allRanks)
        geneIndex <- which(sortRanks == geneRank)
      }
      # add the result to the permutation list
      permGI <- c(permGI, geneIndex)
    }
    # finally add to output: the tool, the patient ID, and the median of the permutation results
    out <- rbind(out, data.frame(tool = df, lovd = row$lovd, gindex = sort(permGI)[3])) # 3rd index for 5 permutations
  }
}

# must be 8 * 308 = 2464
dim(out)[1] == 2464

# do cross-tool count of how many 1-ranked, 2-ranked, etc
outuniq <- as.data.frame(table(out$gindex, out$tool))
colnames(outuniq) <- c("gindex", "tool", "freq")
outuniq$gindex <- as.numeric(as.character(outuniq$gindex))

# add a 0 freq at position 20 to make the plot consistent
for(df in dfs){ outuniq <- rbind(outuniq, data.frame(gindex=20, tool=df, freq=0)) }
  
# fix names ...
levels(outuniq$tool)[match("gado",levels(outuniq$tool))] <- "GADO"
levels(outuniq$tool)[match("vibe",levels(outuniq$tool))] <- "VIBE"
levels(outuniq$tool)[match("amelie",levels(outuniq$tool))] <- "AMELIE"
levels(outuniq$tool)[match("phenomizer",levels(outuniq$tool))] <- "Phenomizer"
levels(outuniq$tool)[match("phenotips",levels(outuniq$tool))] <- "Phenotips"
levels(outuniq$tool)[match("phenix",levels(outuniq$tool))] <- "PhenIX"
levels(outuniq$tool)[match("hiphive",levels(outuniq$tool))] <- "hiPHIVE"
levels(outuniq$tool)[match("pubcf",levels(outuniq$tool))] <- "PubCaseF."

# in addition, do cumulative count of solved cased by rank
outuniq <- outuniq %>%
  group_by(tool) %>%
  mutate(cumuFreq = cumsum(freq))

# visualize and save
ggplot() +
  geom_point(data = outuniq, aes(x = gindex, y = cumuFreq, color = tool), size=3) +
  geom_line(data = outuniq, aes(x = gindex, y = cumuFreq, color = tool), size=1) +
  scale_color_manual(values = colours) +
  scale_x_continuous(breaks = seq(1,20,1), limits = c(1,20)) +
  theme(text = element_text(size=20), legend.title=element_blank(), legend.position = c(0.7, 0.6)) +
  labs(x = "Gene rank in simulated spiked-in clinical gene sets", y = "Number of causal genes detected") +
ggSaveCustom("Figure3", width=8, height=5)


