#!/user/bin/env python3

from shutil import which
from argparse import ArgumentParser
from os.path import isfile
from os.path import isdir
from os import listdir


def readPhenotypes(hpoObo):
    """
    Reads the phenotypes from an .obo file.
    :param hpoObo: the path to the .obo file
    :return: dict with phenotype names as keys and their id as value
    """

    # Strings to identify parts of file.
    termString = '[Term]'
    idString = 'id: '
    nameString = 'name: '
    synonymString = 'synonym: "'

    # Default values for processing.
    phenotypes = {}
    hpoId = None
    name = None

    # Goes through the .obo file.
    for line in open(hpoObo):
        # Resets id and name for new phenotype.
        if line.startswith(termString):
            hpoId = None
            name = None
        # Sets id/name when found.
        elif line.startswith(idString):
            hpoId = line.lstrip(idString).strip()
        elif line.startswith(nameString):
            name = line.lstrip(nameString).strip()
        # Sets a synonym as alternative for name when found on a line.
        elif line.startswith(synonymString):
            name = line.lstrip(synonymString).split('"', 1)[0].strip()

        # If a combination of an id and a name/synonym is stored, saves it to the dictionary.
        # Afterwards, resets name to None so that it won't be triggered by every line (unless a new synonym is found).
        if hpoId is not None and name is not None:
            phenotypes[name] = hpoId
            name = None

    return phenotypes

def retrieveLovdPhenotypes(benchmarkData):
    """
    Retrieves the LOVDs with the phenotypes belonging to each LOVD.
    :param benchmarkData: the file from which the retrieve the LOVDs and their phenotypes
    :return: dict with as keys the LOVDs and a list with phenotypes as value for each key
    """
    # Stores all LOVDs
    lovdPhenotypes = {}

    # Goes through the benchmarking file.
    for i, line in enumerate(open(benchmarkData)):
        # Skips first line (header).
        if i == 0:
            continue

        # Splits the values on their separator.
        line = line.rstrip().split('\t')

        # Checks whether this LOVD was already processed (skips if this is the case).
        if line[0] in lovdPhenotypes.keys():
            continue

        # Adds the phenotypes with their corresponding key.
        lovdPhenotypes[line[0]] = line[4].split(';')

    return lovdPhenotypes


def convertPhenotypeNamesToIds(lovdPhenotypes, phenotypeIdsByName):
    """
    Converts phenotype names into phenotype IDs for a dict containing as values a list of phenotype names.
    :param lovdPhenotypes: the dict for which the names should be replaced with HPO IDs (in the value lists)
    :param phenotypeIdsByName: dict with phenotype names as keys and their id as value
    :return: dict with as keys the LOVDs and a list with HPO IDs as value for each key
    """
    for lovd in lovdPhenotypes.keys():
        phenotypeNames = lovdPhenotypes.get(lovd)
        phenotypeIds = []
        for name in phenotypeNames:
            phenotypeIds.append(phenotypeIdsByName.get(name))
        lovdPhenotypes[lovd] = phenotypeIds

    return lovdPhenotypes


def retrieveAllGenes(hgncFile):
    """
    Retrieves all HGNC gene symbols from the complete HGNC dataset as downloadable from https://www.genenames.org/cgi-bin/statistics
    :param hgncFile: path to HGNC file
    :return: a set with all unique HGNC symbols
    """
    hgncs = set()

    for i, line in enumerate(open(hgncFile)):
        if i > 0:
            hgnc = line.split("\t")[1]
            hgncs.add(hgnc)

    return hgncs


def chunkList(list, chunkSize):
    """
    Cuts a list into chunks and returns these chunks within a list.
    :param list: the list to be chunked
    :param chunkSize: the size of each chunk.
    :return: a list with lists of the defined chunkSize
    see also: https://stackoverflow.com/a/312464
    """
    return [list[i:i + chunkSize] for i in range(0, len(list), chunkSize)]

def validateProgram(commandLine):
    """
    Checks whether the command line (which should contain an application name) refers to a valid path. If it is an
    existing application, does nothing. If the path is not valid, exits with an OSError indicating there is a missing
    application.
    :param commandLine: the name of the application as how it would be called using the command line
    :return:
    """
    if not which(commandLine):
        exit(OSError("Error: " + commandLine + " is not available on this system."))


def fileMergerParser():
    """
    Processes the command line arguments for a parser that only requires an input directory which contains all the files
    and an output file to which the merged output should be written to.
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


def mergeFiles(fileProcessor, inDir, outFile):
    """
    Actual processing of a Python script that merges the files from a single directory into a single file.
    :param fileProcessor: a function that processes a single file stored in the inDir
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

            # Processes file content.
            fileProcessor(fileWriter, filePath)

    # Flushes and closes writer.
    fileWriter.flush()
    fileWriter.close()
