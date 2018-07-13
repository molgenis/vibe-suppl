#!/user/bin/env python3
"""
Name:
    GeneNetworkBenchmarkFileGenerator.py
Example:
    GeneNetworkBenchmarkFileGenerator.py input/ output.tsv

Description:
    Processes the output from AmelieBenchmarkRunner.py for usage in R plots.
"""

from BenchmarkGenerics import fileMergerParser
from BenchmarkGenerics import mergeFiles


def main():
    # Runs application processes.
    args = fileMergerParser()
    mergeFiles(processGeneNetworkFile, args.inDir, args.out)


def processGeneNetworkFile(fileWriter, filePath):
    """
    Processes a single output file from the gene network benchmark.
    :param fileWriter: the file to write the output to
    :param filePath: the path to the file to be processed
    :return:
    """
    # Toggle for first line after the hashtag lines (lines should not be skipped anymore starting from the second
    # non-hashtag line).
    header = True

    # Stores the genes.
    genes = []

    # Writes the genes in order to file separated by a comma (with a newline at the end).
    for line in open(filePath):
        # Skips lines starting with a hashtag.
        if line.startswith("#"):
            continue

        # First line after hashtag lines is the header (describing columns).
        if header:
            header = False
            continue

        # Collects genes.
        genes.append(line.split("\t")[0])

    # Writes genes to file.
    fileWriter.write(",".join(genes))


if __name__ == '__main__':
    main()
