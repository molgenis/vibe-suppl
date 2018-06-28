#!/user/bin/env python3
"""
Name:
    GeneNetworkApiOutputGenerator.py
Example:
    GeneNetworkApiOutputGenerator.py hp.obo benchmark_data.tsv api_output/

Description:
    Retrieves data from the API https://www.genenetwork.nl/api/v1 and writes the genes with weightedZScore to files
    (1 file per LOVD from the benchmark file).
"""

from os.path import isfile
from os.path import isdir
from argparse import ArgumentParser
from time import time
from requests import post
from requests.exceptions import ConnectionError
from requests.exceptions import ReadTimeout
from requests.packages.urllib3 import disable_warnings
from requests.packages.urllib3.exceptions import InsecureRequestWarning
from requests import HTTPError
from BenchmarkGenerics import readPhenotypes
from BenchmarkGenerics import retrieveLovdPhenotypes
from BenchmarkGenerics import convertPhenotypeNamesToIds
from BenchmarkGenerics import waitTillElapsed


def main():
    # Disables InsecureRequestWarning. See also: https://urllib3.readthedocs.org/en/latest/security.html
    disable_warnings(InsecureRequestWarning)

    # Runs application processes.
    args = parseCommandLine()
    phenotypeIdsByName = readPhenotypes(args.hpo)
    lovdPhenotypes = retrieveLovdPhenotypes(args.tsv)
    lovdPhenotypes = convertPhenotypeNamesToIds(lovdPhenotypes, phenotypeIdsByName)
    retrieveGeneNetworkResults(lovdPhenotypes, args.out)


def parseCommandLine():
    """
    Processes the command line arguments.
    :return: args
    """

    # Defines command line.
    parser = ArgumentParser()
    parser.add_argument("hpo", help="he HPO .obo file containing phenotype id's/names")
    parser.add_argument("tsv", help="the benchmarking .tsv file where the first column is the sample ID and the 5th column 1 or more phenotypes (separated by a ';')")
    parser.add_argument("out", help="the file to write output to")

    # Processes command line.
    args = parser.parse_args()

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

    return args


def retrieveGeneNetworkResults(lovdPhenotypes, outDir):
    """
    Retrieves the results from gene network and writes these to files on a per-LOVD basis. If the output dir already
    contains a file with the LOVD name (<lovd>.tsv), that LOVD is skipped (allowing of continuing the benchmark later on
    if stopped).
    :param lovdPhenotypes: benchmark data with as key the LOVD and as value a list of HPO IDs
    :param outDir: the directory to write the output files to (and used to check whether a benchmark for that LOVD was
    already done)
    :return:
    """
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

        # File to write output to.
        fileWriter = open(outFile, 'w')

        # Writes the header to the file.
        fileWriter.write("gene\tweightedZScore\n")

        # The request uri used to retrieve the results.
        uri = "https://www.genenetwork.nl/api/v1/prioritization/" + ",".join(lovdPhenotypes.get(lovd))

        # Waits till elapsed time exceeds 1 second.
        waitTillElapsed(1, time() - requestTime)

        # Tries to make a request to the REST API with the JSON String.
        # If an HTTPError is triggered, this is printed and then no further benchmarking data will be uploaded.
        try:
            response = post(uri, verify=False, timeout=(6, 12))
            response.raise_for_status()
        except (ConnectionError, HTTPError, ReadTimeout) as e:
            exit(e)

        # Stores the current time for managing time between requests.
        requestTime = time()

        # Writes results to file.
        for result in response.json()["results"]:
            fileWriter.write(result["gene"]["name"] + "\t" + str(result["weightedZScore"]) + "\n")

        # Closes file.
        fileWriter.flush()
        fileWriter.close()


if __name__ == '__main__':
    main()
