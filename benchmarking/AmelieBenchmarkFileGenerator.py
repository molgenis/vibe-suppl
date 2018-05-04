#!/user/bin/env python3
"""
Name:
    AmelieBenchmarkFileGenerator.py
Example:
    AmelieBenchmarkFileGenerator.py input/ output.tsv

Description:
    Processes the output from AmelieApiOutputGenerator.py for usage in R plots.
"""

from os.path import isfile
from os.path import isdir
from os import listdir
from argparse import ArgumentParser


def main():
    # Runs application processes.
    args = parseCommandLine()
    mergeFiles(args.inDir, args.out)


def parseCommandLine():
    """
    Processes the command line arguments.
    :return: args
    """

    # Defines command line.
    parser = ArgumentParser()
    parser.add_argument("inDir", help="the directory containing all input files")
    parser.add_argument("out", help="the file to write output to")

    # Processes command line.
    args = parser.parse_args()

    # Removes last slash from directories for consistency.
    args.inDir = args.inDir.rstrip("/")

    # Validates command line.
    if not isdir(args.inDir):
        parser.error('"' + args.inDir.split('/')[-1] + '" is not a valid directory')

    if not args.out.endswith(".tsv"):
        parser.error('"' + args.out.split('/')[-1] + '" is not a .tsv file')
    if isfile(args.out):
        parser.error('"' + args.out.split('/')[-1] + '" already exists')

    return args


def mergeFiles(inDir, outFile):
    """
    Actual processing of the AmelieApiOutputGenerator.py output to a single file with only LOVDs with the ordered genes.
    :param inDir: the directory containing all input files
    :param outFile: the file to write output to
    :return:
    """
    # File to write output to.
    fileWriter = open(outFile, 'w')

    # Writes the header to the file.
    fileWriter.write("lovd\tsuggested_genes\n")

    # Processes all input files.
    for file in listdir(inDir):
        # Full path to file.
        filePath = inDir + "/" + file

        # Skips file if it does not end with ".tsv" or if it is the output file.
        if filePath.endswith(".tsv") and filePath != outFile:
            print("processing: " + file)

            # Writes LOVD to file as first column.
            fileWriter.write(file.split(".")[0] + "\t")

            # Writes the genes in order to file separated by a comma (with a newline at the end).
            for i,line in enumerate(open(filePath)):
                # Skips header line.
                if i == 0:
                    continue
                if i > 1:
                    fileWriter.write(",")
                fileWriter.write(line.split("\t")[0])
            fileWriter.write("\n")

    # Flushes and closes writer.
    fileWriter.flush()
    fileWriter.close()


if __name__ == '__main__':
    main()
