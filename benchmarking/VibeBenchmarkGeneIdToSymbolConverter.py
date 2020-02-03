#!/user/bin/env python3
"""
Name:
    VibeBenchmarkGeneIdToSymbolConverter.py

Example:
    VibeBenchmarkFileGenerator.py fileToProcess.tsv geneInfoFile.csv

Description:
    Converts Gene IDs to symbols for vibe benchmark file (using TSV file generated through DisGeNET as gene ID/symbol source).
"""

from argparse import ArgumentParser
from os.path import isfile

def main():
    # Runs application processes.
    args = fileMergerParser()
    geneInfoDict = storeGeneInfo(args.geneInfoFile)
    convertIdsToSymbols(args.fileToProcess, geneInfoDict)

def fileMergerParser():
    """
    Processes the command line arguments for a parser that only requires an input directory which contains all the files
    and an output file to which the merged output should be written to.
    :return: args
    """

    # Defines command line.
    parser = ArgumentParser()
    parser.add_argument("fileToProcess", help="the file to write output to")
    parser.add_argument("geneInfoFile", help="the directory containing all input files")

    # Processes command line.
    args = parser.parse_args()

    # Validates command line.
    if not isfile(args.fileToProcess):
        parser.error('"' + args.fileToProcess.split('/')[-1] + '" does not exist')
    if not args.fileToProcess.endswith(".tsv"):
        parser.error('"' + args.fileToProcess.split('/')[-1] + '" is not a .tsv file')

    if not isfile(args.geneInfoFile):
        parser.error('"' + args.geneInfoFile.split('/')[-1] + '" does not exist')
    if not args.geneInfoFile.endswith(".csv"):
        parser.error('"' + args.geneInfoFile.split('/')[-1] + '" is not a .csv file')

    return args


def storeGeneInfo(geneInfoFile):
    # Stores the genes with as key the ID and as value the symbol.
    geneInfoDict = {}

    # Goes through all ID/symbol pairs.
    for counter, line in enumerate(open(geneInfoFile)):
        if counter > 0:
            line = line.rstrip().split(",")
            geneInfoDict[line[0]] = line[1]

    return geneInfoDict


def convertIdsToSymbols(fileToProcess, geneInfoDict):
    # File to write output to.
    fileWriter = open(fileToProcess.rpartition(".")[0] + "_converted.tsv", 'w')

    # Goes through all IDs.
    for counter, line in enumerate(open(fileToProcess)):
        # Simply copies the header to new file.
        if counter == 0:
            fileWriter.write(line)
        else:
            # Splits line on separator.
            line = line.rstrip().split("\t")

            # Retrieves the gene IDs.
            geneIds = line[1].split(",")

            # Goes through gene IDs and stores their symbols.
            geneSymbols = []
            for geneId in geneIds:
                geneSymbols.append(geneInfoDict.get(geneId))

            # Writes adjusted line to output file.
            fileWriter.write(line[0] + "\t" + ",".join(geneSymbols) + "\n")

    # Flushes and closes writer.
    fileWriter.flush()
    fileWriter.close()



if __name__ == '__main__':
    main()
