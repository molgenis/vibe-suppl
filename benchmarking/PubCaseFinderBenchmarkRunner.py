#!/user/bin/env python3
"""
Name:
    PubCaseFinderBenchmarkRunner.py
Example:
    PubCaseFinderBenchmarkRunner.py hp.obo benchmark_data.tsv output/

Description:
    Retrieves data from PubCaseFinder using the API public API, processes this
    and writes the output to multiple files (1 file per LOVD from the benchmark file).
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
from BenchmarkGenerics import readPhenotypes
from BenchmarkGenerics import retrieveLovdPhenotypes
from BenchmarkGenerics import convertPhenotypeNamesToIds


def main():
    # Disables InsecureRequestWarning. See also: https://urllib3.readthedocs.org/en/latest/security.html
    disable_warnings(InsecureRequestWarning)

    # Runs application processes.
    args = parseCommandLine()
    phenotypeIdsByName = readPhenotypes(args.hpo)
    lovdPhenotypes = retrieveLovdPhenotypes(args.tsv)
    lovdPhenotypes = convertPhenotypeNamesToIds(lovdPhenotypes, phenotypeIdsByName)
    retrieveExomiserResults(lovdPhenotypes, args.out)


def parseCommandLine():
    """
    Processes the command line arguments.
    :return: args
    """

    # Defines command line.
    parser = ArgumentParser()
    parser.add_argument("hpo", help="The HPO .obo file containing phenotype id's/names")
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


def retrieveExomiserResults(lovdPhenotypes, outDir):
    """
    Retrieves the results from PubCaseFinder and writes these to files on a per-LOVD basis. If the output dir already contains
    a file with the LOVD name (<lovd>.tsv), that LOVD is skipped (allowing of continuing the benchmark later on if stopped).
    :param lovdPhenotypes: benchmark data with as key the LOVD and as value a list of HPO IDs
    :param outDir: the directory to write the output files to (and used to check whether a benchmark for that LOVD was
    already done)
    :return:
    """
    # Goes through all LOVDs.
    for lovd in lovdPhenotypes.keys():
        # Generate phenotype dict list.
        inputPhenotypes = []
        for phenotype in lovdPhenotypes.get(lovd):
            inputPhenotypes.append({"id":phenotype})

        # Defines output file for this LOVD.
        outFile = outDir + "/" + lovd + ".tsv"

        # Checks if output folder already contains a file for this LOVD, and if so, skips this LOVD.
        if isfile(outFile):
            print("# skipping: " + lovd)
            continue

        print("# processing: " + lovd)

        # File to write output to.
        fileWriter = open(outFile, 'w')

        # Tries to make a request to the REST API with the JSON String.
        # If an HTTPError is triggered, this is printed and then no further benchmarking data will be uploaded.
        try:

            response = post("https://pubcasefinder.dbcls.jp/mme/match", verify=False, timeout=(6,600),
                            json={"patient":{"features":inputPhenotypes}})
            response.raise_for_status()
        except (ConnectionError, HTTPError, ReadTimeout) as e:
            exit(e)

        # Digests the results.
        firstGene = True
        for candidate in response.json()['results']:
            genes = candidate['patient']['genomicFeatures']
            for gene in genes:
                if not firstGene:
                    fileWriter.write(',')
                firstGene = False
                fileWriter.write(gene['gene']['label'])

        # Flushes and closes writer.
        fileWriter.flush()
        fileWriter.close()


if __name__ == '__main__':
    main()
