#!/user/bin/env python3
"""
Name:
    BiobesuBenchmarkDataConverter.py
Example:
    BiobesuBenchmarkDataConverter.py hp.obo benchmark_data.tsv output.tsv

Description:
    Converts the benchmark_data.tsv to a "biobesu hpogenerank" compatible format.
"""

from os.path import isfile
from argparse import ArgumentParser
from BenchmarkGenerics import readPhenotypes


def main():
    # Runs application processes.
    args = parseCommandLine()
    phenotypeIdsByName = readPhenotypes(args.hpo)
    convertBenchmarkFile(args.tsv, args.out, phenotypeIdsByName)


def parseCommandLine():
    """
    Processes the command line arguments.
    :return: args
    """

    # Defines command line.
    parser = ArgumentParser()
    parser.add_argument("hpo", help="he HPO .obo file containing phenotype id's/names")
    parser.add_argument("tsv", help="the benchmarking .tsv file where the first column is the sample ID and the 5th column 1 or more phenotypes (separated by a ';')")
    parser.add_argument("out", help="the file to write output to")

    # Processes command line.
    args = parser.parse_args()

    # Validates command line.
    if not args.hpo.endswith(".obo"):
        parser.error('"' + args.hpo.split('/')[-1] + '" is not an .obo file')
    if not isfile(args.hpo):
        parser.error('"' + args.hpo.split('/')[-1] + '" is not an existing file')

    if not args.tsv.endswith(".tsv"):
        parser.error('"' + args.tsv.split('/')[-1] + '" is not a .tsv file')
    if not isfile(args.tsv):
        parser.error('"' + args.tsv.split('/')[-1] + '" is not an existing file')

    if not args.out.endswith(".tsv"):
        parser.error('"' + args.out.split('/')[-1] + '" is not a .tsv file')

    return args

def convertBenchmarkFile(benchmarkData, outFile, phenotypeIdsByName):
    fileWriter = open(outFile, "w")

    # Goes through the benchmarking file.
    for i, line in enumerate(open(benchmarkData)):
        # Skips first line (header).
        if i == 0:
            fileWriter.write("lovd\tgene\thpo_terms\n")
            continue

        # Splits the values on their separator.
        line = line.rstrip().split('\t')

        # Goes through the HPO names and retrieves their code.
        hpoCodes = []
        for hpoName in line[4].split(';'):
            hpoCodes.append(phenotypeIdsByName.get(hpoName))

        # Replaces the HPO names for their codes in the line variable.
        hpoColumn = ';'.join(hpoCodes)

        # Joins line together again.
        fileWriter.write('\t'.join([line[0], line[1], hpoColumn]) + "\n")

    fileWriter.flush()
    fileWriter.close()


if __name__ == '__main__':
    main()
