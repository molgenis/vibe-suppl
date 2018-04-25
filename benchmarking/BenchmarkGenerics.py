#!/user/bin/env python3

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