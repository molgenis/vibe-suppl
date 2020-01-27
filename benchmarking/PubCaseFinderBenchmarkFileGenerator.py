#!/user/bin/env python3
"""
Name:
    PubCaseFinderBenchmarkFileGenerator.py
Example:
    PubCaseFinderBenchmarkFileGenerator.py input/ output.tsv

Description:
    Processes the output from PubCaseFinderBenchmarkRunner.py for usage in R plots.
"""

from BenchmarkGenerics import fileMergerParser
from BenchmarkGenerics import mergeFiles


def main():
    # Runs application processes.
    args = fileMergerParser()
    mergeFiles(processAmelieFile, args.inDir, args.out)


def processAmelieFile(fileWriter, filePath):
    """
    Processes a single output file from the amelie benchmark.
    :param fileWriter: the file to write the output to
    :param filePath: the path to the file to be processed
    :return:
    """
    # Writes the genes in order to file separated by a comma (with a newline at the end).
    fileWriter.write(open(filePath).readline())


if __name__ == '__main__':
    main()
