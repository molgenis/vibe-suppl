#!/usr/bin/env bash

# Sets output file (tmp & final).
OUT_FILE_TMP="$1"/merged.tmp
OUT_FILE_FINAL="$1"/merged.tsv

# printf: prints format of header for file
# >: writes this header to the output file
#
# grep: retrieves full line (each file contains a single line with genes belonging to the LOVD in the filename)
# --- /path/to/<lovd>.tsv:<genes>
# sed: removes full path
# --- <lovd>.tsv:<genes>
# sed: converts ".tsv:" to "\t"
# >>: adds the LOVD with their genes to the output file (format: "LOVD\tgenes" where the genes are separated using commas)
#
# mv: renames the output file (OUT_FILE_TMP -> OUT_FILE_FINAL so that output file would not interfere with the collection
#     of results from the input tsv files)
printf 'lovd\tsuggested_genes\n' > ${OUT_FILE_TMP} && grep '' "$1"/*.tsv | sed "s#^$1/##" | sed $'s/.tsv:/\t/' >> ${OUT_FILE_TMP} && mv ${OUT_FILE_TMP} ${OUT_FILE_FINAL}
