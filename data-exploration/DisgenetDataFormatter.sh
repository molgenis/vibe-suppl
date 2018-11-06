#!/usr/bin/env bash
# Defines error echo.
function errcho { echo -e "$@" 1>&2; }

# Describes usage.
USAGE="Usage: disgenet-data-formatter.sh -p <path> [-h]

Description: Generates tsv files with the unique gene-disease associations from DisGeNET based on \"all_gene_disease_pmid_associations.tsv.gz\" (also includes the non-pmid associations).

Arguments:
-p --path		Path to directory containing the \"all_gene_disease_pmid_associations.tsv.gz\" file.
-h --help Shows 	this help message.
"

function main {
	echo "Digesting command line..."
	digestCommandLine $@
	echo "Generating file containing all associations..."
	digestFile $OUTPUT_FILE_ALL $".*"
	echo "Generating file containing only disease associations..."
	digestFile $OUTPUT_FILE_DISEASE $"^disease"
	echo "Generating file containing only phenotype associations..."
	digestFile $OUTPUT_FILE_PHENOTYPE $"^phenotype$"
	echo "Done."
}

function digestCommandLine {
	# Checks if any arguments were given.
	if [ $# -eq 0 ]; then errcho "No arguments were given.\n\n$USAGE"; exit 1; fi

	#Digests the command line arguments.
	while [[ $# -gt 0 ]]
	do
		key="$1"
		case $key in
			-p|--path)
			DIR="$2"
			shift # argument
			shift # value
			;;
			-h|--help)
			HELP=TRUE
			shift # argument
			;;
			*)    # unknown option
			shift # argument
			;;
		esac
	done

	# Checks if usage is requested.
	if [[ $HELP == TRUE ]]; then echo -e "$USAGE"; exit 0; fi

	# Checks if DIR variable is set. -> http://wiki.bash-hackers.org/syntax/pe#use_an_alternate_value
	if [[ ! ${DIR+isset} == isset ]]; then errcho "Missing required argument: -p/--path <path>\n\n$USAGE"; exit 1; fi

	# Checks if given argument is an existing directory.
	if [ ! -d "$DIR" ]; then errcho "Path is not an existing directory.\n\n$USAGE"; exit 1; fi

	# Checks if directory contains required input file.
	INPUT_FILE=$DIR/all_gene_disease_pmid_associations.tsv.gz
	if [ ! -f "$INPUT_FILE" ]; then errcho "The file \"all_gene_disease_pmid_associations.tsv.gz\" could not be found in the given directory.\n\n$USAGE"; exit 1; fi

	# Sets output files.
	OUTPUT_FILE_ALL=$DIR/complete_gene_associations.tsv
	OUTPUT_FILE_DISEASE=$DIR/complete_gene_disease_associations.tsv
	OUTPUT_FILE_PHENOTYPE=$DIR/complete_gene_phenotype_associations.tsv
}

function digestFile {
	OUTPUT_FILE=$1
	DISEASE_TYPE_REGEX="$2"

	# printf: prints format of header for file
	# >: writes this header to the output file
	#
	# gzip: reads input file
	# --- columns: 1=geneId, 2=diseaseId, 7=originalSource, 9=diseaseType
	# awk: starting from second line (NR>1), print columns [1, 2, 7] if column 9 adheres to regex
	# sort: sort lines
	# uniq: count how many each unique line occurs (adds numbers at the start of the line with trailing spaces)
	# sed: replaces spaces that trailed the uniq results with tabs
	# --- columns: 1=uniq -c output, 2=diseaseId, 3=diseaseId, 4=originalSource
	# awk: merges each unique gene-disease-source line into a gene-disease line
	# --- columns: 1=geneId, 2=diseaseId, 3=total number of associations, 4=number of sources, 5=associations per source (format: <source1>=<count1>;<source2>=<count2>;<...>)
	# sort: sorts lines again
	# >>: adds output to file containing the header
	#
	# gzip: gzips the file (if gzipped file already exists, it is overwritten)
	printf 'geneId\tdiseaseId\tnAssociations\tnSources\tsources\n' > $OUTPUT_FILE && gzip -cd $INPUT_FILE | awk -v diseaseTypeRegex="$DISEASE_TYPE_REGEX" -F $'\t' 'NR>1 {if ($9 ~ diseaseTypeRegex) print $1 "\t" $2 "\t" $7}' | sort | uniq -c | sed -e 's/ *//' -e $'s/ /\t/' | awk -F $'\t' '{nAssociations[$2"\t"$3]+=$1; nSources[$2"\t"$3]++; sources[$2"\t"$3]=sources[$2"\t"$3] ";" $4 "=" $1} END{for(geneDisease in nAssociations) print geneDisease "\t" nAssociations[geneDisease] "\t" nSources[geneDisease] "\t" substr(sources[geneDisease],2)}' | sort >> $OUTPUT_FILE && gzip -f $OUTPUT_FILE
}

main $@
