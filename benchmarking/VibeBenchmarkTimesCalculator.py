#!/user/bin/env python3
"""
Name:   VibeBenchmarTimesCalculator.py

Example:
    VibeBenchmarTimesCalculator.py errorLogDir/ outputFile.tsv

Description:
    TODO.
"""

from argparse import ArgumentParser
from os.path import isfile
from os.path import isdir
from os import listdir
from re import match
from math import floor

def main():
    # Runs application processes.
    args = parseCommandLine()
    calculateTimes(args.logDir, args.out)


def parseCommandLine():
    """
    Processes the command line arguments
    :return: args
    """

    # Defines command line.
    parser = ArgumentParser()
    parser.add_argument("logDir", help="directory containing the error logs (which contain information produced by the BSD time tool)")
    parser.add_argument("out", help="the file to write output to")

    # Processes command line.
    args = parser.parse_args()

    # Removes last slash from directories for consistency.
    args.logDir = args.logDir.rstrip("/")

    # Validates command line.
    if not isdir(args.logDir):
        parser.error('"' + args.logDir.split('/')[-1] + '" is not a valid directory')

    if not args.out.endswith(".tsv"):
        parser.error('"' + args.out.split('/')[-1] + '" is not a .tsv file')
    if isfile(args.out):
        parser.error('"' + args.out.split('/')[-1] + '" already exists')

    return args

def calculateTimes(inDir, outFile):
    # File to write output to.
    fileWriter = open(outFile, 'w')

    # Writes the header to the file.
    fileWriter.write("id\ttime (mm:ss)\n")

    # Processes all input files.
    for file in sorted(listdir(inDir)):
        # Full path to file.
        inFilePath = inDir + "/" + file

        if inFilePath.endswith(".log"):
            print("processing: " + inFilePath)

            # Writes identifier to file as first column.
            fileWriter.write(file.split(".")[0].lstrip("err_") + "\t")

            # Stores total time (user+sys) for file.
            timeNeeded = 0.0

            # Process lines in file.
            for line in open(inFilePath):
                # If user or sys is found, retrieves time and stores it in timeNeeded as seconds (float).
                if line.startswith("user") or line.startswith("sys"):
                    line = line.split("\t")[1]
                    lineMatch = match("(\d+)m(\d+(.\d+)?)s", line)
                    timeNeeded += float(lineMatch.group(1)) * 60 + float(lineMatch.group(2))

            # After user and sys times were collected, converts seconds back to minutes and (whole) seconds.
            minutes = floor(timeNeeded / 60)
            seconds = round(timeNeeded % 60)

            # Writes time to file (m:s).
            fileWriter.write("{:02d}:{:02d}".format(minutes, seconds))

            # Writes a newline so that each processed file has its own line in the output file.
            fileWriter.write("\n")

    # Flushes and closes writer.
    fileWriter.flush()
    fileWriter.close()

if __name__ == '__main__':
    main()