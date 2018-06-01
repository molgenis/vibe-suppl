# vibe-suppl
Variant Interpretation using Biological Evidence - supplemental files

## Info
This repo contains supplemental files regarding the Java application found [here][vibe].

## Benchmarking

There are several benchmarking scripts available with some generic code used by multiple benchmarks in a separate file.
An explanation on how to run the can be found below.

* __`AmelieApiOutputGenerator.py`__
    * __Info:__ Connects to `https://amelie.stanford.edu/api/` to retrieve the gene scores for each set of HPO terms
    available in the benchmark data. As the genes of interest should be entered manually and there is a limit in the
    number of entered genes, the [complete HGNC dataset][hgnc_complete]
    is used and divided over multiple separate requests so that all genes get a score. As the scores are only sorted
    per request, a sort on all genes is done prior to file writing.
    * __Example:__ `AmelieApiOutputGenerator.py hp.obo hgnc_complete_set.txt benchmark_data.tsv output/`
* __`AmelieBenchmarkFileGenerator.py`__
    * __Info:__ Converts the output from `AmelieApiOutputGenerator.py` for usage in `BenchmarkResultsProcessor.R`.
    * __Example:__ `AmelieBenchmarkFileGenerator.py input/ output.tsv`
* __`BenchmarkGenerics.py`__
    * __Info:__ Contains methods used in multiple scripts.
* __`BenchmarkResultsProcessor.R`__
    * __Info:__ Creates plots from the benchmark data.
* __`PhenomizerBenchmarkRunner.py`__
    * __Info:__ Uses the [query_phenomizer][query_phenomizer] python tool to process all benchmark data.
    * __Important:__ Before running this tool, be sure to install `query_phenomizer` using the following code:
    ```
    git clone https://github.com/svandenhoek/query_phenomizer.git
    cd query_phenomizer
    pip install --editable .
    ``` 
    * __Example:__ `PhenomizerBenchmarkRunner.py username hp.obo benchmark_data.tsv output/`
* __`PhenotipsBenchmarkFileGenerator.py`__
    * __Info:__ Uses the API of Phenotips to upload the benchmark dataset and then download the results.
    * __Important:__ A phenotips instince to which can be connected is required. Please refer to the
    [Phenotips download page][phenotips_download] for more information.
    * __Example:__ `PhenotipsBenchmarkFileGenerator.py http://localhost:8080/ username hp.obo benchmark_data.tsv out.tsv`
* __`VibeBenchmarkBashScriptsGenerator.py`__
    * __Info:__ Generates bash files used for benchmarking. For each
    * __Important:__ As each VIBE instance needs a separate database, please refer to the information in the script
    itself for how to prepare for the benchmarking correctly.
    * __Example:__ `VibeBenchmarkBashScriptsGenerator.py hp.obo benchmark_data.tsv scripts/ 50`

There are several files used among these scripts. These include:
* benchmark_data.tsv
    * A dataset with the first column being an ID and the fourth column 1 or more phenotypes separated
    by a comma (the phenotype names should exist within the [Human Phenotype Ontology][hpo_obo]) .
* [hp.obo][hpo_obo]
    * The Human Phenotype Ontology used for combining/converting phenotype names with their HPO ID.
* [hgnc_complete_set.txt][hgnc_complete]
    * The HUGO Gene Nomenclature Committee file containing information about genes (primarily used to generate a list
    containing all genes).

[vibe]:https://github.com/molgenis/vibe
[hgnc_complete]:ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
[query_phenomizer]:https://github.com/svandenhoek/query_phenomizer
[phenotips_download]:https://phenotips.org/Download
[hpo_obo]:http://purl.obolibrary.org/obo/hp.obo