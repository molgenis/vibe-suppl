#!/user/bin/env python3
"""
Name:
    AmelieBenchmarkFileGenerator.py
Example:
    AmelieBenchmarkFileGenerator.py api_output/ output.tsv

Description:
    Processes the output from AmelieApiOutputGenerator.py for usage in R plots.
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
