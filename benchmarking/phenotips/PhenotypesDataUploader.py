#!/user/bin/env python3
"""
Name: PhenotypesDataUploader.py
Usage: PhenotypesDataUploader.py <benchmark TSV file> <HPO OBO file>

Example:
    PhenotypesDataUploader.py benchmark_data.tsv hp.obo


Description:

"""

from sys import argv
from requests import post
from requests import HTTPError
from requests.auth import HTTPBasicAuth


def main():
    phenotypes = readPhenotypes(argv[2])
    uploadPhenotypes(argv[1], phenotypes)


def readPhenotypes(hpoObo):
    termString = '[Term]'
    idString = 'id: '
    nameString = 'name: '
    synonymString = 'synonym: "'

    phenotypes = {}
    hpoId = None
    name = None

    for line in open(hpoObo):
        if line.startswith(termString):
            hpoId = None
            name = None
        elif line.startswith(idString):
            hpoId = line.lstrip(idString)
        elif line.startswith(nameString):
            name = line.lstrip(nameString)
        elif line.startswith(synonymString):
            name = line.lstrip(synonymString).split('"', 1)[0]

        if hpoId is not None and name is not None:
            phenotypes[name.strip()] = hpoId.strip()
            name = None

    return phenotypes


def uploadPhenotypes(dataToUpload, phenotypes):
    for i, line in enumerate(open(dataToUpload)):
        if i == 0:
            continue

        line = line.rstrip().split('\t')

        requestString = '{"solved":{"status":"unsolved"},"external_id":"' + line[0] + '","clinicalStatus":"affected","features":['

        for j, hpoName in enumerate(line[4].split(';')):
            if j > 0:
                requestString += ','

            requestString += '{"id":"' + phenotypes.get(hpoName) + '","label":"' + hpoName + '","type":"phenotype","observed":"yes"}'

        requestString += "]}"

        try:
            response = post('http://localhost:8080/rest/patients', data=requestString, auth=HTTPBasicAuth('Admin', 'admin'))
            response.raise_for_status()
        except HTTPError as e:
            print(e)


if __name__ == '__main__':
    main()
