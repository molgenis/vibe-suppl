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

from BenchmarkFileGeneSymbolToIdConverter import parseCommandLine
from BenchmarkFileGeneSymbolToIdConverter import storeGeneInfo
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

            # Goes through gene IDs and stores their symbols.
            possibleNewValues = geneInfoDict.get(line[1].upper())
            if possibleNewValues is None:
                eprint("No replacement NCBI gene id could be found for {}, using symbol instead.".format(line[1]))
            elif possibleNewValues[0] != '':
                line[1] = possibleNewValues[0]
            elif possibleNewValues[1] != '':
                #print("Using uncurated match for {}".format(oldValue))
                line[1] = possibleNewValues[1]
            else:
                raise Exception("An unknown exception has occured for {}.".format(line[1]))

            # Writes adjusted line to output file.
            fileWriter.write("\t".join(line) + "\n")

    # Flushes and closes writer.
    fileWriter.flush()
    fileWriter.close()



if __name__ == '__main__':
    main()
