#!/user/bin/env python3
"""
Name:
    BenchmarkFileGeneSymbolToIdConverter.py

Example:
    BenchmarkFileGeneSymbolToIdConverter.py fileToProcess.tsv geneInfoFile.csv

Description:
    Converts HGNC gene symbols to NCBI gene IDs. Uses a file generated using:
    https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_prev_sym&col=md_eg_id&col=gd_pub_eg_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_hgnc_id&format=text&submit=submit

    Expected header:
    HGNC ID	Approved symbol	Previous symbols	NCBI Gene ID(supplied by NCBI)	NCBI Gene ID
"""

from argparse import ArgumentParser
from os.path import isfile
from os import linesep
from sys import stderr

def eprint(*args, **kwargs):
    """
    Error printing function.
    :param args:
    :param kwargs:
    :return:
    """
    print(*args, file=stderr, **kwargs)

def main():
    # Runs application processes.
    args = parseCommandLine()
    geneInfoDict = storeGeneInfo(args.geneInfoFile)
    convertIdsToSymbols(args.fileToProcess, geneInfoDict)

def parseCommandLine():
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
    if not args.geneInfoFile.endswith(".tsv"):
        parser.error('"' + args.geneInfoFile.split('/')[-1] + '" is not a .tsv file')

    return args


def storeGeneInfo(geneInfoFile):
    # Expected header to validate if input file is correct.
    expectedHeader = "HGNC ID\tApproved symbol\tPrevious symbols\tNCBI Gene ID(supplied by NCBI)\tNCBI Gene ID"
    # Stores the genes with as key the symbol and as value an array with curated gene ID and NCBI provided gene ID.
    geneInfoDictApproved = {}
    geneInfoDictPrevious = {}

    # Goes through all approved gene symbols.
    for counter, line in enumerate(open(geneInfoFile)):
        # Validates if all expected columns are present and in expected order.
        if counter == 0 and line.rstrip() != expectedHeader:
            raise Exception("Unexpected gene info file header.{0}Expected: {1}{0}Actual: {2}".format(linesep, expectedHeader, line))

        # Processes items.
        if counter > 0:
            line = line.rstrip().split("\t")

            # Ensures all items are of equal length after splitting.
            while len(line) < 5:
                line.append('')

            # Validates if curated and NCBI provided genes are equal if both are provided.
            if line[4] != '' and line[3] != '' and line[4] != line[3]:
                eprint("Non-matching gene IDs for {}: {} according to HGNC and {} according to NCBI.".format(line[1], line[4], line[3]))

            # Adds approvedSymbol:(curatedId,providedId) to dict if either an approved or provided gene ID is available.
            if line[4] != '' or line[3] != '':
                if geneInfoDictApproved.get(line[1]) is not None:
                    raise Exception("{} was already mentioned as approved symbol with provided gene IDs.".format(line[1]))
                else:
                    addItemToGeneInfoDict(geneInfoDictApproved, line[1], line)

    # Goes through all previous gene symbols.
    for counter, line in enumerate(open(geneInfoFile)):
        # Processes items.
        if counter > 0:
            line = line.rstrip().split("\t")

            # Ensures all items are of equal length after splitting.
            while len(line) < 5:
                line.append('')

            if line[4] != '' or line[3] != '':
                for prevSymbol in line[2].split(","):
                    if prevSymbol != "":
                        if geneInfoDictApproved.get(prevSymbol) is not None:
                            print("Gene IDs were already found for {} as approved symbol. Skipping it as previous symbol.".format(prevSymbol))
                        elif geneInfoDictPrevious.get(prevSymbol) is not None:
                            eprint("Gene IDs were already found for {} as previous symbol. Skipping additional previous symbol result.".format(prevSymbol))
                        else:
                            addItemToGeneInfoDict(geneInfoDictPrevious, prevSymbol.strip(), line)

    # Merges approved with previous symbols (approved overrides previous, though this should not happen anyhow).
    geneInfoDictPrevious.update(geneInfoDictApproved)
    return geneInfoDictPrevious

def addItemToGeneInfoDict(geneInfoDict, symbol, line):
    geneInfoDict[symbol.upper()] = (line[4], line[3])

def convertIdsToSymbols(fileToProcess, geneInfoDict):
    # File to write output to.
    fileWriter = open(fileToProcess.rpartition(".")[0] + "_ncbi-gene-ids.tsv", 'w')

    # Goes through all IDs.
    for counter, line in enumerate(open(fileToProcess)):
        # Simply copies the header to new file.
        if counter == 0:
            fileWriter.write(line)
        else:
            # Splits line on separator.
            line = line.rstrip().split("\t")

            # Simply writes the line if no results were generated.
            if len(line) == 1:
                fileWriter.write(line[0] + "\t" + linesep)
                continue

            # Retrieves the gene IDs.
            oldValues = line[1].split(",")

            # Goes through gene IDs and stores their symbols.
            newValues = []
            for oldValue in oldValues:
                possibleNewValues = geneInfoDict.get(oldValue.upper())
                if possibleNewValues is None:
                    eprint("No replacement NCBI gene id could be found for {}, using old symbol instead.".format(oldValue))
                    newValues.append(oldValue)
                elif possibleNewValues[0] != '':
                    newValues.append(possibleNewValues[0])
                elif possibleNewValues[1] != '':
                    #print("Using uncurated match for {}".format(oldValue))
                    newValues.append(possibleNewValues[0])
                else:
                    raise Exception("An unknown exception has occured for {}.".format(oldValues))


            # Writes adjusted line to output file.
            fileWriter.write(line[0] + "\t" + ",".join(newValues) + "\n")

    # Flushes and closes writer.
    fileWriter.flush()
    fileWriter.close()



if __name__ == '__main__':
    main()
