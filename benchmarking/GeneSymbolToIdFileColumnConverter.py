#!/user/bin/env python3
"""
Name:
    GeneSymbolToIdFileColumnConverter.py

Example:
    GeneSymbolToIdFileColumnConverter.py fileToProcess.tsv geneInfoFile.csv

Description:
    Converts HGNC gene symbols to NCBI gene IDs. Uses a file generated using:
    https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_prev_sym&col=md_eg_id&col=gd_pub_eg_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_hgnc_id&format=text&submit=submit

    Expected header:
    HGNC ID	Approved symbol	Previous symbols	NCBI Gene ID(supplied by NCBI)	NCBI Gene ID
"""

from BenchmarkFileGeneSymbolToIdConverter import storeGeneInfo
from argparse import ArgumentParser
from os.path import isfile
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
    convertSymbolsToIds(args.fileToProcess, geneInfoDict, args.column, args.skipNoMatch)

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
    parser.add_argument("column", type=int, help="the column of the file that stores the HGNC gene symbol")
    parser.add_argument("skipNoMatch", type=int, help="whether a result where no match was found should skipped in the output (0: uses symbol instead, 1:skips line)")

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

    if args.skipNoMatch not in [0,1]:
        parser.error('skipNoMatch can only be 0 or 1')

    return args

def convertSymbolsToIds(fileToProcess, geneInfoDict, columnNumber, skipNoMatch):
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
            possibleNewValues = geneInfoDict.get(line[columnNumber].upper())
            if possibleNewValues is None:
                if skipNoMatch:
                    eprint("No replacement NCBI gene id could be found for {} (line {}), skipping line.".format(line[columnNumber], counter))
                    continue
                else:
                    eprint("No replacement NCBI gene id could be found for {} (line {}), using symbol instead.".format(line[columnNumber], counter))
            elif possibleNewValues[0] != '':
                line[columnNumber] = possibleNewValues[0]
            elif possibleNewValues[1] != '':
                #print("Using uncurated match for {}".format(oldValue))
                line[columnNumber] = possibleNewValues[1]
            else:
                raise Exception("An unknown exception has occured for {}.".format(line[columnNumber]))

            # Writes adjusted line to output file.
            fileWriter.write("\t".join(line) + "\n")

    # Flushes and closes writer.
    fileWriter.flush()
    fileWriter.close()



if __name__ == '__main__':
    main()
