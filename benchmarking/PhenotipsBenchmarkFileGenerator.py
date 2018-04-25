#!/user/bin/env python3
"""
Name: PhenotipsBenchmarkFileGenerator.py

Example:
    PhenotipsBenchmarkFileGenerator.py http://localhost:8080/ Admin hp.obo benchmark_data.tsv out.tsv

Description:
    Uploads benchmark data to Phenotips using HPO data (as an .obo file) and then downloads suggested genes.
"""

from os.path import isfile
from argparse import ArgumentParser
from requests import post
from requests import get
from requests.exceptions import ConnectionError
from requests import HTTPError
from requests.auth import HTTPBasicAuth
from getpass import getpass
from re import match
from BenchmarkGenerics import readPhenotypes
from BenchmarkGenerics import retrieveLovdPhenotypes


def main():
    args = parseCommandLine()
    phenotypeIdsByName = readPhenotypes(args.hpo)
    lovdPhenotypes = retrieveLovdPhenotypes(args.tsv)
    pwd = getpass("Phenotips password for " + args.username + ":")
    uploadPhenotypes(args.url, args.username, pwd, phenotypeIdsByName, lovdPhenotypes)
    downloadGenes(args.url, args.username, pwd, args.out, lovdPhenotypes.keys())


def parseCommandLine():
    """
    Processes the command line arguments
    :return: args
    """

    # Defines command line.
    parser = ArgumentParser()
    parser.add_argument("url", help="the url to the phenotips instance to upload to (")
    parser.add_argument("username", help="the username for authentication to the phenotips server")
    parser.add_argument("hpo", help="he HPO .obo file containing phenotype id's/names")
    parser.add_argument("tsv", help="the benchmarking .tsv file where the first column is the sample ID and the 5th column 1 or more phenotypes (separated by a ';')")
    parser.add_argument("out", help="the file to write output to")

    # Processes command line.
    args = parser.parse_args()

    # Strips slash on the end of an url.
    args.url = args.url.rstrip("/")

    # Validates command line.
    if match("https?:\/\/[a-z0-9.\-]+(:[0-9]{4})?", args.url) is None:
        parser.error("invalid url")
    if match("^\w+$", args.username) is None:
        parser.error("the username may only contain: greek characters (a-z A-Z), numbers (0-9) and underscores (_)")
    if not args.hpo.endswith(".obo"):
        parser.error('"' + args.hpo.split('/')[-1] + '" is not an .obo file')
    if not isfile(args.hpo):
        parser.error('"' + args.hpo.split('/')[-1] + '" is not an existing file')
    if not args.tsv.endswith(".tsv"):
        parser.error('"' + args.tsv.split('/')[-1] + '" is not a .tsv file')
    if not isfile(args.tsv):
        parser.error('"' + args.tsv.split('/')[-1] + '" is not an existing file')
    if not args.out.endswith(".tsv"):
        parser.error('"' + args.out.split('/')[-1] + '" is not a .tsv file')
    if isfile(args.out):
        parser.error('"' + args.out.split('/')[-1] + '" already exists')

    return args


def uploadPhenotypes(phenotipsUrl, username, password, phenotypeIdsByName, lovdPhenotypes):
    """
    Uploads the benchmark data to phenotips.
    :param phenotipsUrl: the url to upload the benchmark data to
    :param username: the username for authentication
    :param password: the password for authentication
    :param phenotypeIdsByName: dict with phenotype names as keys and their id as value
    :param lovdPhenotypes: benchmark data with as key the LOVD and as value a list of phenotype names
    :return: list with all LOVDs that were uploaded
    """

    for lovd in lovdPhenotypes.keys():
        # Starts generating the JSON for the request.
        requestString = '{"solved":{"status":"unsolved"},"external_id":"' + lovd + '","clinicalStatus":"affected","features":['

        # Goes through all phenotypes and adds these to the JSON String.
        for phenotypeName in lovdPhenotypes.get(lovd):
            requestString += '{"id":"' + phenotypeIdsByName.get(phenotypeName) + '","label":"' + phenotypeName + '","type":"phenotype","observed":"yes"}'
        requestString += "]}"

        # Tries to make a request to the REST API with the JSON String.
        # If an HTTPError is triggered, this is printed and then no further benchmarking data will be uploaded.
        try:
            response = post(phenotipsUrl + "/rest/patients", data=requestString, auth=HTTPBasicAuth(username, password))
            response.raise_for_status()
        except (ConnectionError, HTTPError) as e:
            exit(e)


def downloadGenes(phenotipsUrl, username, password, out, lovds):
    """
    Downloads the suggested genes for each uploaded LOVD
    :param phenotipsUrl: the url to upload the benchmark data to
    :param username: the username for authentication
    :param password: the password for authentication
    :param out: the file to write the results to
    :param lovds: list with the LOVDs for which the suggested genes should be retrieved
    """

    # File to write output to.
    fileWriter = open(out, 'w')

    # Writes the header to the file.
    fileWriter.write("lovd\tsuggested_genes\n")

    # Goes through all LOVDs.
    for lovd in lovds:
        # Writes the LOVD followed by a tab
        fileWriter.write(lovd + "\t")

        # Tries to make a several request to the REST API for data retrieval.
        # If an HTTPError is triggered, this is printed and then no further benchmarking data will be uploaded.
        try:
            # Retrieves internal ID based on external ID.
            response = get(phenotipsUrl + "/rest/patients/eid/" + lovd, auth=HTTPBasicAuth(username, password))
            response.raise_for_status()

            # If external ID is only used once, directly retrieve internal ID.
            try:
                phenotipsId = response.json()['report_id']
            # Otherwise retrieve the internal ID from the first item.
            # Any items that have the same external ID are assumed to be equal.
            except KeyError:
                phenotipsId = response.json()['patients'][0]['id']

            # Retrieves the suggested genes for the LOVD.
            response = get(phenotipsUrl + "/rest/patients/" + phenotipsId + "/suggested-gene-panels", auth=HTTPBasicAuth(username, password))
            response.raise_for_status()
        except (ConnectionError, HTTPError) as e:
            exit(e)

        # Goes through all suggested genes for a single LOVD.
        for i, gene in enumerate(response.json()['genes']):
            # If more than 1 gene found, adds separator to output between genes.
            if i > 0:
                fileWriter.write(',')
            fileWriter.write((gene['gene_symbol']))

        # After all suggested genes are processed, adds a newLine.
        fileWriter.write('\n')

    # Flushes and closes writer.
    fileWriter.flush()
    fileWriter.close()


if __name__ == '__main__':
    main()
