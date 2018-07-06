#!/user/bin/env python3
"""
Name:   VibeBenchmarkBashScriptGenerator.py

Example:
    VibeBenchmarkBashScriptsGenerator.py hp.obo benchmark_data.tsv scripts/ 50

Description:
    Generates bash scripts for running the vibe benchmark. Assumes that for each bash script a unique TDB is present
    (this is because a TDB can only be accessed by 1 JVM at a time). The scripts should be ran in a folder with the
    following structure:

    + main_folder
    |- vibe_benchmark_0.sh
    |- vibe_benchmark_1.sh
    |- ...
    |
    |- TDB_0/
    |- TDB_1/
    |- ...
    |
    |- vibe-with-dependencies.jar
    |
    |- results/
    |- out/
    |- err/

    Where "..." indicates more files like the one from the line above.

    The number given to the application indicates the maximum number of runs done by a single bash script
     (50 in the example).
"""

from argparse import ArgumentParser
from os.path import isfile
from os.path import isdir
from BenchmarkGenerics import readPhenotypes
from BenchmarkGenerics import retrieveLovdPhenotypes


def main():
    args = parseCommandLine()
    phenotypeIdsByName = readPhenotypes(args.hpo)
    lovdPhenotypes = retrieveLovdPhenotypes(args.tsv)
    generateBashScripts(phenotypeIdsByName, lovdPhenotypes, args.out, args.runLimit)


def parseCommandLine():
    """
    Processes the command line arguments
    :return: args
    """

    # Defines command line.
    parser = ArgumentParser()
    parser.add_argument("hpo", help="he HPO .obo file containing phenotype id's/names")
    parser.add_argument("tsv", help="the benchmarking .tsv file where the first column is the sample ID and the 5th column 1 or more phenotypes (separated by a ';')")
    parser.add_argument("out", help="the directory to write output to")
    parser.add_argument("runLimit", help="the number of benchmark java calls per bash file")

    # Processes command line.
    args = parser.parse_args()

    # Removes last slash from directories for consistency.
    args.out = args.out.rstrip("/")

    # Converts input values to correct type (throws error for invalid type).
    try:
        args.runLimit = int(args.runLimit)
    except ValueError:
        parser.error('runLimit must be a number')

    # Validates command line.
    if not args.hpo.endswith(".obo"):
        parser.error('"' + args.hpo.split('/')[-1] + '" is not an .obo file')
    if not isfile(args.hpo):
        parser.error('"' + args.hpo.split('/')[-1] + '" is not an existing file')
    if not args.tsv.endswith(".tsv"):
        parser.error('"' + args.tsv.split('/')[-1] + '" is not a .tsv file')
    if not isfile(args.tsv):
        parser.error('"' + args.tsv.split('/')[-1] + '" is not an existing file')
    if not isdir(args.out):
        parser.error('"' + args.out.split('/')[-1] + '" is not a valid directory')
    if args.runLimit < 1:
        parser.error('runLimit must be 1 or higher: ' + args.runLimit)

    return args


def generateBashScripts(phenotypeIdsByName, lovdPhenotypes, outDir, runLimit):
    """
    Generates the bash files.
    :param phenotypeIdsByName: dict with phenotype names as keys and their id as value
    :param lovdPhenotypes: dict with as keys the LOVDs and a list with phenotypes as value for each key
    :param outDir: the directory to write the bash files to
    :param runLimit: the limit in number subsequent java calls per bash script
    :return:
    """
    fileCounter = 0  # Keeps track of the number of files/filenames.
    runCounter = 0  # Keeps track of the number of subsequent java calls for a single file.

    # Creates initial bash file.
    bashFile = createNewWritableBashFile(outDir, fileCounter)

    # Goes through all LOVDs
    for lovd in sorted(lovdPhenotypes.keys()):
        # If limit is reached: finishes file, updates counters and generates a new file.
        if runCounter >= runLimit:
            bashFile.flush()
            bashFile.close()
            fileCounter += 1
            runCounter = 0
            bashFile = createNewWritableBashFile(outDir, fileCounter)

        # If not first java call in bash script, adds code for between calls.
        if runCounter > 0:
            bashFile.write(" && \n")

        # Writes the basic part of a java call.
        bashFile.write("time( java -XX:ParallelGCThreads=1  -Xmx4g -jar vibe-with-dependencies.jar -v -t TDB_" + str(fileCounter) + "/ -s gda_max -o results/" + lovd + ".tsv")

        # Adds all phenotypes to the java call belonging to the LOVD.
        for phenotypeName in lovdPhenotypes.get(lovd):
            bashFile.write(" -p " + phenotypeIdsByName.get(phenotypeName))

        # Finishing part of a java call.
        bashFile.write(" ) 1>out/out_" + lovd + ".log 2>err/err_" + lovd + ".log")
        runCounter += 1

    # Closes last file to which data was written.
    bashFile.flush()
    bashFile.close()


def createNewWritableBashFile(outDir, fileCounter):
    """
    Generates a new bash file to write commands to.
    :param outDir: the directory to write the bash files to
    :param fileCounter: number describing which file it it
    :return:
    """
    # Creates file.
    bashFile = open(outDir + "/vibe_benchmark_" + str(fileCounter) + ".sh", "w")

    # Writes initial line.
    bashFile.write("#!/usr/bin/env bash\n")

    # Returns file for writing.
    return bashFile


if __name__ == '__main__':
    main()