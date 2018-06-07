# BenchmarkingResultsProcessor.R

########
# Name:
# singleBenchmarkResultProcessor
#
# Description:
# Digests a single benchmark
#
# Input:
# benchmarkDataRow -  A single row from the benchmark dataset (used for running
#                     the benchmarks). The first column should contain an LOVD
#                     while the second one a gene.
#                     Should be 1-dimensional (f.e. using apply or unlist)!
# resultData -  The data.frame containing the results from a full benchmark.
#               The row names should contain the LOVDs while the first column
#               should be the genes (in the retrieved order) separated by commas.
#
# Output:
# 
########
singleBenchmarkResultProcessor <- function(benchmarkDataRow, resultData) {
  match(benchmarkDataRow[2],
        strsplit(resultData[benchmarkDataRow[1],], ",")[[1]])
}

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

phenotips <- read.table(paste0(baseDir,"results/phenotips.tsv"), header=T,
                        sep="\t", colClasses=c("character"), row.names=1)
phenomizer <- read.table(paste0(baseDir,"results/phenomizer.tsv"),
                         header=T, sep="\t", colClasses=c("character"),
                         row.names=1)
amelie <- read.table(paste0(baseDir,"results/amelie.tsv"),
                     header=T, sep="\t", colClasses=c("character"),
                     row.names=1)
vibe.20180503 <- read.table(paste0(baseDir,"results/vibe_2018-05-02.tsv"),
                            header=T, sep="\t", colClasses=c("character"),
                            row.names=1)

# Process data
phenotips.results <-apply(benchmarkData, 1, singleBenchmarkResultProcessor,
                          resultData=phenotips)

phenomizer.results <-apply(benchmarkData, 1, singleBenchmarkResultProcessor,
                              resultData=phenomizer)

amelie.results <-apply(benchmarkData, 1, singleBenchmarkResultProcessor,
                           resultData=amelie)

vibe.20180503.results <-apply(benchmarkData, 1, singleBenchmarkResultProcessor,
                              resultData=vibe.20180503)


allResults <- data.frame(amelie.results, phenomizer.results, phenotips.results,
                         vibe.20180503.results)
colnames(allResults) <- c("amelie", "phenomizer", "phenotips", "vibe")

# Count missing values.
resultsWithAnNa <- allResults[apply(apply(allResults, c(1,2), is.na), 1, any),]
sum(apply(apply(resultsWithAnNa, c(1,2), is.na), 1, all))
naCounts <- apply(apply(resultsWithAnNa, c(1,2), is.na), 2, sum)

# Create plot
yAxisMax <- ceiling(max(allResults, na.rm=T)/100)*100

postscript(paste0(imgExportDir, 'benchmarking_comparison.eps'), width=6, height=8)
boxplot(allResults, xaxt='n',yaxt='n', pch="")
abline(h=1, col="gray92")
abline(h=seq(100, yAxisMax, 100), col="gray92")
boxplot(allResults, las=1, pch=20, yaxt='n',
        main="position of relevant genes among different tools",
        xlab="tool", ylab="position", col="white", add=T)
axis(2, las=1, at=1)
axis(2, las=1, at=seq(500, yAxisMax,500))
dev.off()

