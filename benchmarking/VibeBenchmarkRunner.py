#!/user/bin/env python3
"""
Name: VibeBenchmarkRunner.py

Example:
    VibeBenchmarkRunner.py vibe-with-dependencies.jar ./TDB/ hp.owl hp.obo benchmark_data.tsv 0 ./output/

Description:

"""

from argparse import ArgumentParser
from os.path import isfile
from os.path import isdir
from BenchmarkGenerics import readPhenotypes
from BenchmarkGenerics import retrieveLovdPhenotypes
from copy import copy
import subprocess


def main():
    args = parseCommandLine()
    phenotypeIdsByName = readPhenotypes(args.hpoObo)
    lovdPhenotypes = retrieveLovdPhenotypes(args.tsv)
    runBenchmark(phenotypeIdsByName, lovdPhenotypes, args.basicJarArgs, args.out)


def parseCommandLine():
    """
    Processes the command line arguments
    :return: args
    """

    # Defines command line.
    parser = ArgumentParser()
    parser.add_argument("jar", help="the vibe application .jar file")
    parser.add_argument("tdb", help="the DisGeNET TDB")
    parser.add_argument("hpoOwl", help="the HPO .owl file containing phenotype id's/names")
    parser.add_argument("hpoObo", help="he HPO .obo file containing phenotype id's/names")
    parser.add_argument("tsv", help="the benchmarking .tsv file where the first column is the sample ID and the 5th column 1 or more phenotypes (separated by a ';')")
    parser.add_argument("maxDistance", help="the max distance value to use (if 0, app will be run without phenotypes retrieval phase)")
    parser.add_argument("out", help="the directory to write output to")

    # Processes command line.
    args = parser.parse_args()

    # Converts input values to correct type (throws error for invalid type).
    try:
        args.maxDistance = int(args.maxDistance)
    except ValueError:
        parser.error('maxDistance must be a number')

    # Removes last slash from directories for consistency.
    args.tdb = args.tdb.rstrip("/")
    args.out = args.out.rstrip("/")

    # Validates command line.
    if not args.jar.endswith(".jar"):
        parser.error('"' + args.jar.split('/')[-1] + '" is not a .jar file')
    if not isfile(args.jar):
        parser.error('"' + args.jar.split('/')[-1] + '" is not an existing file')
    # Currently tdb validation itself handled by jar
    if not isdir(args.tdb):
        parser.error('"' + args.tdb.split('/')[-1] + '" is not a valid directory')
    if not args.hpoOwl.endswith(".owl"):
        parser.error('"' + args.hpoOwl.split('/')[-1] + '" is not an .owl file')
    if not isfile(args.hpoOwl):
        parser.error('"' + args.hpoOwl.split('/')[-1] + '" is not an existing file')
    if not args.hpoObo.endswith(".obo"):
        parser.error('"' + args.hpoObo.split('/')[-1] + '" is not an .obo file')
    if not isfile(args.hpoObo):
        parser.error('"' + args.hpoObo.split('/')[-1] + '" is not an existing file')
    if not args.tsv.endswith(".tsv"):
        parser.error('"' + args.tsv.split('/')[-1] + '" is not a .tsv file')
    if args.maxDistance < 0:
        parser.error('maxDistance must be 0 or higher: ' + args.maxDistance)
    if not isfile(args.tsv):
        parser.error('"' + args.tsv.split('/')[-1] + '" is not an existing file')
    if not isdir(args.out):
        parser.error('"' + args.out.split('/')[-1] + '" is not a valid directory')

    # Processes maxDistance.
    if args.maxDistance == 0:
        args.basicJarArgs = ["java", "-jar", args.jar, "-v", "-t", args.tdb]
    else:
        args.basicJarArgs = ["java", "-jar", args.jar, "-v", "-t", args.tdb, "-w", args.hpoOwl, "-d", "-m", str(args.maxDistance)]

    # Validates the output directory whether it was used with a different benchmark configuration and if no config.txt
    # file is present creates one for this run.
    if not benchmarkConfigFileProcessor(args.basicJarArgs, args.out):
        parser.error('the given output directory was used with a different benchmark configuration')


    return args


def benchmarkConfigFileProcessor(basicJarArgs, outputDir):
    """
    Validates the output directory whether it was used with a different benchmark configuration and if no config.txt
    file is present creates one with the configuration of this benchmark run.
    :param basicJarArgs: the basic java arguments to run a jar and arguments specific for the benchmark run
    :param outputDir: the directory to write the output files to
    :return: True if non-used dir or uses same config, false if a different config was detected
    """
    # The config file.
    argsFile = outputDir + "/" + "config.txt"

    # A String version of the basicJarArgs.
    basicJarArgsString = " ".join(basicJarArgs)

    # If config file found checks content.
    if isfile(argsFile):
        argsFileReader = open(argsFile)
        line = argsFileReader.readline()
        argsFileReader.close()
        if line != basicJarArgsString:
            return False
    # If no config file found, creates one.
    else:
        argsFileWriter = open(argsFile, "w")
        argsFileWriter.write(basicJarArgsString)
        argsFileWriter.flush()
        argsFileWriter.close()

    return True


def runBenchmark(phenotypeIdsByName, lovdPhenotypes, basicJarArgs, outputDir):
    """
    Runs the benchmark. If an output file already exists in the output dir, skips this benchmark (in case benchmark was
    abborted and continued later on)
    :param phenotypeIdsByName: dict with phenotype names as keys and their id as value
    :param lovdPhenotypes: dict with as keys the LOVDs and a list with phenotypes as value for each key
    :param basicJarArgs: the basic java arguments to run a jar and arguments specific for the benchmark run
    :param outputDir: the directory to write the output files to
    :return:
    """
    for lovd in sorted(lovdPhenotypes.keys()):
        jarArgs, outFile = prepareArguments(phenotypeIdsByName, basicJarArgs, lovd, lovdPhenotypes.get(lovd), outputDir)
        print(outFile)
        # Skips a run if output file already found (allows for continuing an unfinished benchmark run)
        if isfile(outFile):
            print("########## Skipping benchmark: " + lovd)
        else:
            print("########## Running benchmark: " + lovd)
            runProcess(jarArgs)


def prepareArguments(phenotypeIdsByName, basicJarArgs, lovd, phenotypes, outputDir):
    """
    Prepares the subprocess arguments for all LOVDs
    :param phenotypeIdsByName: dict with phenotype names as keys and their id as value
    :param basicJarArgs: the basic java arguments to run a jar and arguments specific for the benchmark run
    :param lovd: a single LOVD name
    :param phenotypes: the phenotypes belonging to the LOVD specified in parameter "name"
    :param outputDir: the directory to write the output files to
    :return: jarArgs, outFile
                jarArgs: the full command line to be passed to subprocess
                outFile: the output file which the subprocess will write output to
    """
    # Makes a copy of the arguments so that adjustments won't effect original list.
    jarArgs = copy(basicJarArgs)

    # Creates and adds output file to arguments.
    outFile = outputDir + "/" + lovd + ".tsv"
    jarArgs += ["-o", outFile]

    # Adds all phenotypes to the arguments.
    for phenotypeName in phenotypes:
        jarArgs += ["-p", phenotypeIdsByName.get(phenotypeName)]

    return jarArgs, outFile


def runProcess(args):
    """
    Runs a single subprocess and continues to write stdout until the subprocess has stopped.
    :param jarArgs: the arguments used for subprocess
    :return:
    """
    # Starts a process.
    process = subprocess.Popen(args)

    # Keeps looking for stdout from the process until it stopped.
    while process.poll() is None:
        if process.stdout is not None:
            print(process.stdout.readline().rstrip())


if __name__ == '__main__':
    main()
