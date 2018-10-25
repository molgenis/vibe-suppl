#!/user/bin/env python3
"""
Name:
    ToolsBenchmarkGeneScoresCollector.py
Example:
    ToolsBenchmarkGeneScoresCollector.py benchmarkdata.tsv outfile.tsv amelie_output/ geneNetwork_output/ vibe_output/

Description:
    Retrieves highest score per gene from the different tool benchmarks.
"""

from argparse import ArgumentParser
from os.path import isfile
from os.path import isdir
from re import match

def main():
    args = commandLineParser()
    genesForIds = retrieveBenchmarkIds(args.benchmarkdata)
    generateFileWithScores(args, genesForIds)

def commandLineParser():
    """
    Processes the command line arguments for a parser that only requires an input directory which contains all the files
    and an output file to which the merged output should be written to.
    :return: args
    """

    # Defines command line.
    parser = ArgumentParser()
    parser.add_argument("benchmarkdata", help="the benchmark data file used")
    parser.add_argument("outfile", help="the file to write output to")
    parser.add_argument("amelie", help="the directory with amelie input files")
    parser.add_argument("genenetwork", help="the directory with gene network input files")
    #parser.add_argument("phenomizer", help="the directory with phenomizer input files")
    #parser.add_argument("phenotips", help="the directory with phenotips input files")
    parser.add_argument("vibe", help="the directory with vibe input files")

    # Processes command line.
    args = parser.parse_args()

    # Removes last slash from directories for consistency.
    args.amelie = args.amelie.rstrip("/")
    args.genenetwork = args.genenetwork.rstrip("/")
    args.vibe = args.vibe.rstrip("/")

    # Validates command line.
    if not isfile(args.benchmarkdata):
        parser.error('"' + args.benchmarkdata.split('/')[-1] + '" does not exist')
    if not args.benchmarkdata.endswith(".tsv"):
        parser.error('"' + args.benchmarkdata.split('/')[-1] + '" is not a .tsv file')

    if not args.outfile.endswith(".tsv"):
        parser.error('"' + args.outfile.split('/')[-1] + '" is not a .tsv file')
    if isfile(args.outfile):
        parser.error('"' + args.outfile.split('/')[-1] + '" already exists')

    if not isdir(args.amelie):
        parser.error('"' + args.amelie.split('/')[-1] + '" is not a valid directory')
    if not isdir(args.genenetwork):
        parser.error('"' + args.genenetwork.split('/')[-1] + '" is not a valid directory')
    if not isdir(args.vibe):
        parser.error('"' + args.vibe.split('/')[-1] + '" is not a valid directory')

    return args

def retrieveBenchmarkIds(benchmarkData):
    """
    Retrieves patient/case IDs with the validation gene(s) belonging to that patient/case.
    :param benchmarkData: the benchmark data file
    :return: a dict with as key the patient/case ID and as value a list with genes belonging to it
    """
    genesForIds = {}

    for i, line in enumerate(open(benchmarkData)):
        # Skips first line (header).
        if i == 0:
            continue

        lineSplits = line.split("\t")
        id = lineSplits[0]
        gene = lineSplits[1]

        if id in genesForIds.keys():
            genesForIds.get(id).append(gene)
        else:
            genesForIds[id] = [gene]

    return genesForIds


def generateFileWithScores(args, genesForIds):
    """
    Generates a file with the scores for each patient/case-gene combination.
    :param args: the command line arguments
    :param genesForIds: a dict with as key the patient/case ID and as value a list with genes belonging to it
    :return:
    """
    # File to write output to.
    fileWriter = open(args.outfile, 'w')

    # Writes the header to the file.
    fileWriter.write("lovd\tgene\tamelie\tgeneNetwork\tvibe\n")

    # Goes through all the patient/case IDs and for each of these goes through all genes.
    for id in genesForIds.keys():
        print("Processing ID: " + id)
        for gene in genesForIds.get(id):
            # Retrieves score for gene (or stores the specified value if gene was not found in output).
            amelieScore = retrieveScoreFromTool(id, gene, args.amelie, "^([a-zA-Z0-9\-]+)\t", "^[a-zA-Z0-9\-]+\t[0-9]+:([0-9.]+)", "0")
            geneNetworkScore = retrieveScoreFromTool(id, gene, args.genenetwork, "^([a-zA-Z0-9\-]+)\t", "^[a-zA-Z0-9\-]+\t[a-zA-Z0-9]+\t[0-9]+\t([0-9.\-]+)\t", "NA")
            #retrieveScoreFromTool(id, gene, args.phenomizer)
            #retrieveScoreFromTool(id, gene, args.phenotips)
            vibeScore = retrieveScoreFromTool(id, gene, args.vibe, "^([a-zA-Z0-9\-]+)\t", "^[a-zA-Z0-9\-]+\t.*\t([0-9.]+)\t[0-9.]+\t[0-9.]+\n", "0")

            # Writes the tool scores with this patient/case-gene combination to the file.
            fileWriter.write(id + "\t" + gene + "\t" + amelieScore + "\t" + geneNetworkScore + "\t" + vibeScore + "\n")

    fileWriter.flush()
    fileWriter.close()

def retrieveScoreFromTool(id, gene, toolDir, geneFieldRegex, scoreFieldRegex, naValue):
    """
    Looks for the gene in the tool output and returns it if possible. If not, naValue is returned instead.
    :param id: the patient/case ID (used as fileName in the output files from the tools)
    :param gene: the gene to look for
    :param toolDir: the directory containing the output for this tool (which stores the patient/case ID files)
    :param geneFieldRegex: the regular expression used to find if the line represents the gene which is looked for
    :param scoreFieldRegex: the regular expression used to extract the score belonging to a gene
    :param naValue: the value to return if this file does not contain the gene which is looked for
    :return: if the gene was found in the file, the score belonging to it. Otherwise, the given naValue
    """
    for i, line in enumerate(open(toolDir + "/" + id + ".tsv")):
        geneFieldRegexMatch = match(geneFieldRegex, line)
        if geneFieldRegexMatch and geneFieldRegexMatch.group(1) == gene:
            return match(scoreFieldRegex, line).group(1)
    return naValue


if __name__ == '__main__':
    main()