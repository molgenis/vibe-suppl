#!/user/bin/env python3
"""
Name:
    AmelieApiOutputGenerator.py
Example:
    AmelieApiOutputGenerator.py hp.obo hgnc_complete_set.txt benchmark_data.tsv output/

Description:
    Retrieves data from https://amelie.stanford.edu/ using the API, processes this and writes the output to multiple
    files (1 file per LOVD from the benchmark file).
"""

from os.path import isfile
from os.path import isdir
from argparse import ArgumentParser
from requests import post
from requests.exceptions import ConnectionError
from requests.exceptions import ReadTimeout
from requests.packages.urllib3 import disable_warnings
from requests.packages.urllib3.exceptions import InsecureRequestWarning
from requests import HTTPError
from time import time
from time import sleep
from BenchmarkGenerics import readPhenotypes
from BenchmarkGenerics import retrieveLovdPhenotypes
from BenchmarkGenerics import convertPhenotypeNamesToIds
from BenchmarkGenerics import retrieveAllGenes
from BenchmarkGenerics import chunkList


def main():
    # Disables InsecureRequestWarning. See also: https://urllib3.readthedocs.org/en/latest/security.html
    disable_warnings(InsecureRequestWarning)

    # Runs application processes.
    args = parseCommandLine()
    phenotypeIdsByName = readPhenotypes(args.hpo)
    lovdPhenotypes = retrieveLovdPhenotypes(args.tsv)
    lovdPhenotypes = convertPhenotypeNamesToIds(lovdPhenotypes, phenotypeIdsByName)
    hgncs = retrieveAllGenes(args.hgnc)
    retrieveAmelieResults(lovdPhenotypes, hgncs, args.out)


def parseCommandLine():
    """
    Processes the command line arguments
    :return: args
    """

    # Defines command line.
    parser = ArgumentParser()
    parser.add_argument("hpo", help="he HPO .obo file containing phenotype id's/names")
    parser.add_argument("hgnc", help="complete HGNC dataset as downloaded from https://www.genenames.org/cgi-bin/statistics")
    parser.add_argument("tsv", help="the benchmarking .tsv file where the first column is the sample ID and the 5th column 1 or more phenotypes (separated by a ';')")
    parser.add_argument("out", help="the file to write output to")

    # Processes command line.
    args = parser.parse_args()

    # Validates command line.
    if not args.hpo.endswith(".obo"):
        parser.error('"' + args.hpo.split('/')[-1] + '" is not an .obo file')
    if not isfile(args.hpo):
        parser.error('"' + args.hpo.split('/')[-1] + '" is not an existing file')

    if not args.hgnc.endswith(".txt"):
        parser.error('"' + args.hgnc.split('/')[-1] + '" is not a .txt file')
    if not isfile(args.hgnc):
        parser.error('"' + args.hgnc.split('/')[-1] + '" is not an existing file')

    if not args.tsv.endswith(".tsv"):
        parser.error('"' + args.tsv.split('/')[-1] + '" is not a .tsv file')
    if not isfile(args.tsv):
        parser.error('"' + args.tsv.split('/')[-1] + '" is not an existing file')

    if not isdir(args.out):
        parser.error('"' + args.out.split('/')[-1] + '" is not a valid directory')

    return args


def retrieveAmelieResults(lovdPhenotypes, hgncs, outDir):
    """
    Retrieves the results from amelie and writes these to files on a per-LOVD basis. If the output dir already contains
    a file with the LOVD name (<lovd>.tsv), that LOVD is skipped (allowing of continuing the benchmark later on if stopped).
    :param lovdPhenotypes: benchmark data with as key the LOVD and as value a list of HPO IDs
    :param hgncs: a set with all unique HGNC symbols
    :return:
    """
    # Chunks for amelie retrieval (too many for all at once).
    hgncs = chunkList(list(hgncs), 1000)

    # Stores initial time as negative time() so that sleep is not triggered the first time.
    requestTime = -time()

    # Goes through all LOVDs.
    for lovd in lovdPhenotypes.keys():
        # Defines output file for this LOVD.
        outFile = outDir + "/" + lovd + ".tsv"

        # Checks if output folder already contains a file for this LOVD, and if so, skips this LOVD.
        if isfile(outFile):
            print("# skipping: " + lovd)
            continue

        print("# processing: " + lovd)

        # Storage of gene results for a single LOVD.
        lovdAmelieOutput = {}

        # Goes through all gene chunks (not all genes can be processed at once).
        for i, hgncsChunk in enumerate(hgncs):
            print("chunk: " + str(i+1) + "/" + str(len(hgncs)))
            # Checks if previous request finished at least 1 second ago, and if not, sleeps 1 second before making
            # another request.
            elapsedTime = time() - requestTime
            if elapsedTime < 1:
                sleep(1)

            # Tries to make a request to the REST API with the JSON String.
            # If an HTTPError is triggered, this is printed and then no further benchmarking data will be uploaded.
            try:
                response = post("https://amelie.stanford.edu/api/", verify=False, timeout=(6,180),
                                data={"genes":",".join(hgncsChunk), "phenotypes":",".join(lovdPhenotypes.get(lovd))})
                response.raise_for_status()
            except (ConnectionError, HTTPError, ReadTimeout) as e:
                exit(e)

            # Stores the current time for managing time between requests.
            requestTime = time()

            # Digests the results for the genes from a single chunk.
            for gene in response.json():
                geneSymbol = gene[0]
                scores = []
                pubmedIds = []

                # Goes through all score/pubmed ID sets.
                for pubmedScores in gene[1]:
                    scores.append(float(pubmedScores[0]))  # Expects float so enforces this.
                    pubmedIds.append(pubmedScores[1])

                # If there were any results for the gene, adds these to the LOVD results.
                if len(scores) > 0:
                    lovdAmelieOutput[geneSymbol] = (scores, pubmedIds)

        # Writes the LOVD output to a file.
        writeLovdResultsToFile(lovdAmelieOutput, outFile)
        exit(0)


def retrieveSortedAmelieList(lovdAmelieOutput):
    """
    Sorts the genes based on their first stored score and returns this as a list. Assumes that the first score is the highest.
    :param lovdAmelieOutput: a dict with as keys a gene symbol and as values a tuple with lists
            (first list contains all scores from highest to lowest, second list contains all pubmed IDs).
    :return: a list with gene symbols sorted based on the first item of the value lists
    """
    return sorted(lovdAmelieOutput, key=lambda k: lovdAmelieOutput[k][0][0], reverse=True)


def writeLovdResultsToFile(lovdAmelieOutput, outFile):
    """
    Writes results for a single LOVD to a file.
    :param lovdAmelieOutput: a dict with as keys a gene symbol and as values a tuple with lists
            (first list contains all scores from highest to lowest, second list contains all pubmed IDs).
    :param outFile: the file to write the output to
    :return:
    """
    # File to write output to.
    fileWriter = open(outFile, 'w')

    # Writes the header to the file.
    fileWriter.write("gene\tscores\tpubmed IDs\n")

    # Goes through all genes and writes these to the file (1 line per gene).
    for gene in retrieveSortedAmelieList(lovdAmelieOutput):
        fileWriter.write(gene + "\t" + ";".join(map(str,lovdAmelieOutput.get(gene)[0])) + "\t" + ";".join(lovdAmelieOutput.get(gene)[1]) + "\n")

    # Flushes and closes writer.
    fileWriter.flush()
    fileWriter.close()


if __name__ == '__main__':
    main()
