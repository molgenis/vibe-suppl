#!/user/bin/env python3
"""
Name:
    BasicBenchmarkFilesMerger.py
Example:
    BasicBenchmarkFilesMerger.py api_output/ output.tsv

Description:
    Processes the output from benchmark scripts which generate output .tsv files containing 1 gene per line with the
    gene being the first column (with optionally containing other columns with additional information).
"""

from BenchmarkGenerics import fileMergerParser
from BenchmarkGenerics import mergeFiles


def main():
    # Runs application processes.
    args = fileMergerParser()
    mergeFiles(processAmelieFile, args.inDir, args.out)


def processAmelieFile(fileWriter, filePath):
    # Writes the genes in order to file separated by a comma (with a newline at the end).
    for i, line in enumerate(open(filePath)):
        # Skips header line.
        if i == 0:
            continue
        if i > 1:
            fileWriter.write(",")
        fileWriter.write(line.split("\t")[0])
    fileWriter.write("\n")


if __name__ == '__main__':
    main()
