#!/user/bin/env python3
"""
Name: PhenotipsBenchmarkFileGenerator.py

Example:
    PhenotipsBenchmarkFileGenerator.py http://localhost:8080/ Admin hp.obo benchmark_data.tsv out.tsv

Description:
    Uploads benchmark data to Phenotips using HPO data (as an .obo file).
"""

from os.path import isfile
from argparse import ArgumentParser
from requests import post
from requests import get
from requests import HTTPError
from requests.auth import HTTPBasicAuth
from getpass import getpass
from re import match


def main():
    args = parseCommandLine()
    phenotypes = readPhenotypes(args.hpo)
    pwd = getpass("Phenotips password for " + args.username + ":")
    lovds = uploadPhenotypes(args.url, args.username, pwd, phenotypes, args.tsv)
    downloadGenes(args.url, args.username, pwd, args.out, lovds)


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
            hpoId = line.lstrip(idString)
        elif line.startswith(nameString):
            name = line.lstrip(nameString)
        # Sets a synonym as alternative for name when found on a line.
        elif line.startswith(synonymString):
            name = line.lstrip(synonymString).split('"', 1)[0]

        # If a combination of an id and a name/synonym is stored, saves it to the dictionary.
        # Afterwards, resets name to None so that it won't be triggered by every line (unless a new synonym is found).
        if hpoId is not None and name is not None:
            phenotypes[name.strip()] = hpoId.strip()
            name = None

    return phenotypes


def uploadPhenotypes(phenotipsUrl, username, password, phenotypes, dataToUpload):
    """
    Uploads the benchmark data to phenotips.
    :param phenotipsUrl: the url to upload the benchmark data to
    :param username: the username for authentication
    :param password: the password for authentication
    :param phenotypes: dict with phenotype names as keys and their id as value
    :param dataToUpload: the benchmark .tsv file
    :return: list with all LOVDs that were uploaded
    """

    # Stores all LOVDs
    lovds = []

    # Goes through the benchmarking file.
    for i, line in enumerate(open(dataToUpload)):
        # Skips first line (header).
        if i == 0:
            continue

        # Splits the values on their separator.
        line = line.rstrip().split('\t')

        # Checks whether this is the same test individual as the previous one and skips it if this is the case.
        # Reason: some individual are stored multiple times as they have multiple gene/OMIM matches.
        if i > 1 and line[0] == lovds[-1]:
            continue

        # Adds the LOVD to the list with all processed LOVDs.
        lovds.append(line[0])

        # Starts generating the JSON for the request.
        requestString = '{"solved":{"status":"unsolved"},"external_id":"' + line[0] + '","clinicalStatus":"affected","features":['

        # Goes through all phenotypes (after being split on their separator) and adds these to the JSON String.
        for j, hpoName in enumerate(line[4].split(';')):
            if j > 0:
                requestString += ','

            requestString += '{"id":"' + phenotypes.get(hpoName) + '","label":"' + hpoName + '","type":"phenotype","observed":"yes"}'

        requestString += "]}"

        # Tries to make a request to the REST API with the JSON String.
        # If an HTTPError is triggered, this is printed and then no further benchmarking data will be uploaded.
        try:
            response = post(phenotipsUrl + "/rest/patients", data=requestString, auth=HTTPBasicAuth(username, password))
            response.raise_for_status()
        except HTTPError as e:
            print(e)
            break

    return lovds


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
        except HTTPError as e:
            print(e)
            break

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
