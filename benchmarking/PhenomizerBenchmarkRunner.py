#!/user/bin/env python3
"""
Name: PhenomizerBenchmarkRunner.py

Example:
    PhenomizerBenchmarkRunner.py username hp.obo benchmark_data.tsv output/

Description:
    Runs the benchmark data through Phenomizer and writes the results to files in the defined output directory.

    Make sure https://github.com/svandenhoek/query_phenomizer is installed (create clone -> cd to directory ->
    "pip install --editable .")!
"""

from argparse import ArgumentParser
from os.path import isfile
from os.path import isdir
from re import match
from re import sub
from getpass import getpass
from BenchmarkGenerics import readPhenotypes
from BenchmarkGenerics import retrieveLovdPhenotypes
from BenchmarkGenerics import convertPhenotypeNamesToIds
from BenchmarkGenerics import validateProgram
from subprocess import call
from time import sleep


def main():
    # Checks if query_phenomizer is available and exits if not.
    validateProgram("query_phenomizer")

    args = parseCommandLine()
    phenotypeIdsByName = readPhenotypes(args.hpo)
    lovdPhenotypes = retrieveLovdPhenotypes(args.tsv)
    lovdPhenotypes = convertPhenotypeNamesToIds(lovdPhenotypes, phenotypeIdsByName)
    pwd = getpass("Phenomizer password for " + args.username + ":")
    retrievePhenomizerResults(lovdPhenotypes, args.username, pwd, args.out)


def parseCommandLine():
    """
    Processes the command line arguments
    :return: args
    """

    # Defines command line.
    parser = ArgumentParser()
    parser.add_argument("username", help="the username for authentication to the phenotips server")
    parser.add_argument("hpo", help="he HPO .obo file containing phenotype id's/names")
    parser.add_argument("tsv", help="the benchmarking .tsv file where the first column is the sample ID and the 5th column 1 or more phenotypes (separated by a ';')")
    parser.add_argument("out", help="the directory to write output to")

    parser.description = "Please make sure that the following code has been run before using this script: " \
                         "git clone https://github.com/svandenhoek/query_phenomizer && cd query_phenomizer && pip install --editable ."

    # Processes command line.
    args = parser.parse_args()

    # Strips slash on the end of an url.
    args.out = args.out.rstrip("/")

    # Validates command line.
    if match("^.+$", args.username) is None:
        parser.error("the username should have at least a length of 1 character.")
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

    return args


def retrievePhenomizerResults(lovdPhenotypes, username, password, outDir):
    """
    Runs the query_phenomizer tool for each LOVD.

    :param lovdPhenotypes: benchmark data with as key the LOVD and as value a list of HPO IDs
    :param username: the username for Phenomizer
    :param password: the password for Phenomizer
    :param outDir: the directory to store the results in
    :return:
    """

    # Goes through all LOVDs.
    for lovd in lovdPhenotypes.keys():
        # Defines output file for this LOVD.
        outFile = outDir + "/" + lovd + ".tsv"

        # Checks if output folder already contains a file for this LOVD, and if so, skips this LOVD.
        if isfile(outFile):
            print("# skipping: " + lovd)
            continue

        print("# processing: " + lovd)

        # Generates the command line to run query_phenomizer.
        commandLine = ["query_phenomizer", "-u", username, "-p", password, "-o", outFile] + lovdPhenotypes.get(lovd)

        # Runs query_phenomizer.
        try:
            call(commandLine)
        except Exception as e:
            processQueryPhenomizerException(e)

        # Sleeps before another run is done so that the server is not overloaded.
        sleep(4)


def processQueryPhenomizerException(exception):
    """
    Prints the error message and exits the application while filtering the password from the output if this output
    contains the line used for running the query_phenomizer.
    :param exception: the exception that is thrown
    :return:
    """
    # Prints exception message. If part of it includes the command line, removes the password so this is not shown.
    print(sub("-p .+ ", "-p <redacted> ", exception.args[1]))
    # Exits the application with the exit code belonging to the exception.
    exit(exception.args[0])


if __name__ == '__main__':
    main()
