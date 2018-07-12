#!/user/bin/env python3
"""
Name:
    PhenomizerBenchmarkFileGenerator.py
Example:
    PhenomizerBenchmarkFileGenerator.py api_output/ output.tsv

Description:
    Processes the output from PhenomizerBenchmarkFileGenerator.py for usage in R plots.
"""

from BenchmarkGenerics import fileMergerParser
from BenchmarkGenerics import mergeFiles
from BenchmarkGenerics import removeDuplicates


def main():
    # Runs application processes.
    args = fileMergerParser()
    mergeFiles(processPhenomizerFile, args.inDir, args.out)


def processPhenomizerFile(fileWriter, filePath):
    # Stores all the genes from the file.
    allGenes = []

    # Writes the genes in order to file separated by a comma (with a newline at the end).
    for i, line in enumerate(open(filePath)):
        # Skips header line.
        if i == 0:
            continue

        # Processes line.
        lineSplits = line.split("\t")

        # If last column is not empty, processes line.
        if len(lineSplits[3].strip()) > 0:
            # Splits the genes and removes the number of the gene in brackets after the name.
            genes = lineSplits[3].split(",")
            for i,gene in enumerate(genes):
                allGenes.append(gene.split("(")[0].strip())

    # Writes the genes to the file while removing duplicates (first mention is kept).
    fileWriter.write(",".join(removeDuplicates(allGenes)))


if __name__ == '__main__':
    main()
