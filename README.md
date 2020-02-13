# vibe-suppl
This repo contains supplemental files regarding the Java application found [here][vibe]. Note that these are in no way
needed to use the vibe tool, but were used to generate additional information (such as benchmarking). They were created
with the assumption that they are used exactly in the way they are meant to be used, so while certain checks/validations
might be present, using these scripts in the wrong way might result in weird behavior.

## Benchmarking

### Scripts

There are several benchmarking scripts available with some generic code used by multiple benchmarks in a separate file.
An explanation on how to run the can be found below. In general, the `Runner` scripts runs the benchmark while the
`FileGenerator` script (if available) formats the `Runner` output to a more usable format. Some exceptions are present,
such as for vibe where there is a `ParallelBashScriptsGenerator` instead. So please refer to to
<a href="#running-the-benchmarks">this section<a/> for more information regarding running the individual benchmarks.

* __`AmelieApiOutputGenerator.py`__
    * __Info:__ Connects to `https://amelie.stanford.edu/api/` to retrieve the gene scores for each set of HPO terms
    available in the benchmark data. As the genes of interest should be entered manually and there is a limit in the
    number of entered genes, the [complete HGNC dataset][hgnc_complete]
    is used and divided over multiple separate requests so that all genes get a score. As the scores are only sorted
    per request, a sort on all genes is done prior to file writing.
* __`AmelieBenchmarkRunner.py`__
    * __Info:__ Converts the output from `AmelieApiOutputGenerator.py` for usage in `BenchmarkResultsProcessor.R`.
* __`BenchmarkFileHpoConverter.py`__
    * __Info:__ A script to convert a benchmark file containing HPO names in the fifth column to a benchmark file with
     HPO codes in the fifth column. Should not be needed for running existing benchmarks, but is supplied as a
     convenience script in case benchmarks are created that cannot use `BenchmarkGenerics.py` but do need HPO codes as
     input. 
* __`BenchmarkGenerics.py`__
    * __Info:__ Contains methods used in multiple scripts.
    * __Important:__ This script should not be ran independently. If Python scripts are moved (for example to a server
    to run the benchmarks there), be sure to include this file within the same directory.
* __`BenchmarkResultsProcessor.R`__
    * __Info:__ Creates plots from the benchmark data.
* __`GeneNetworkBenchmarkFileGenerator.py`__
    * __Info:__ Converts the output from `GeneNetworkBenchmarkRunner.py` for usage in `BenchmarkResultsProcessor.R`.
* __`GeneNetworkBenchmarkRunner.py`__
    * __Info:__ Connects to the API from `https://www.genenetwork.nl/` to retrieve the prioritized genes based on input
    phenotypes.
* __`PhenomizerBenchmarkFileGenerator.py`__
    * __Info:__ Converts the output from `PhenomizerBenchmarkRunner.py` for usage in `BenchmarkResultsProcessor.R`.
* __`PhenomizerBenchmarkRunner.py`__
    * __Info:__ Uses the [query_phenomizer][query_phenomizer] python tool to process all benchmark data.
    * __Important:__ [query_phenomizer][query_phenomizer] needs to be installed on the system. Additionally, an account
    is needed for running [query_phenomizer][query_phenomizer].
* __`PhenotipsBenchmarkRunner.py`__
    * __Info:__ Uses the API of Phenotips to upload the benchmark dataset and then download the results.
    * __Important:__ A phenotips instince to which can be connected is required. Please refer to the
    [Phenotips download page][phenotips_download] for more information.
* __`VibeBenchmarkFileGenerator.py`__
    * __Info:__ Converts the output from `VibeBenchmarkParallelBashScriptsGenerator.py` for usage in `BenchmarkResultsProcessor.R`.
* __`VibeBenchmarkParallelBashScriptsGenerator.py`__
    * __Info:__ Generates bash files used for benchmarking (by using a limit of runs per file). Note that for each
    created bash script a separate TDB is needed. Please refer to the documentation in the script itself for more
    information.
    * __Important:__ As each VIBE instance needs a separate database, please refer to the information in the script
    itself for how to prepare for the benchmarking correctly.

### Data

There are several files used among these scripts. These include:
* benchmark_data.tsv
    * A dataset with the first column being an ID and the fourth column 1 or more phenotypes separated
    by a comma (the phenotype names should exist within the [Human Phenotype Ontology][hpo_obo]) .
* [hp.obo][hpo_obo]
    * The Human Phenotype Ontology used for combining/converting phenotype names with their HPO ID.
* [hgnc_complete_set.txt][hgnc_complete]
    * The HUGO Gene Nomenclature Committee file containing information about genes (primarily used to generate a list
    containing all genes).

### Running the benchmarks

#### Amelie

1. Run benchmark:
    ```
    python3 AmelieBenchmarkRunner.py hp.obo hgnc_complete_set.txt benchmark_data.tsv amelie_output/
    ```

2. Process benchmark output:
    ```
    python3 AmelieBenchmarkFileGenerator.py amelie_output/ amelie_results.tsv
    ```

#### GADO

We used the stand-alone commandline version GADO (v 1.0.1), available at: https://github.com/molgenis/systemsgenetics/wiki/GADO-Command-line. We accepted all automatically suggested alternative HPO terms in cases that the supplied HPO term could not be used. We have used the prediction matrix `hpo_predictions_sigOnly_spiked_01_02_2018`.

#### Phenomizer

1. Install [query_phenomizer][query_phenomizer] (if not already installed):
    ```
    git clone https://github.com/svandenhoek/query_phenomizer.git
    cd query_phenomizer
    pip install --editable .
    ```

2. Run benchmark:
    ```
    python3 PhenomizerBenchmarkRunner.py username hp.obo benchmark_data.tsv phenomizer_output/
    ```

3. Process benchmark output:
    ```
    python3 PhenomizerBenchmarkFileGenerator.py phenomizer_output/ phenomizer_results.tsv
    ```

#### Phenotips

1. Install [phenotips][phenotips_download].

2. Run benchmark:
    ```
    python3 PhenotipsBenchmarkRunner.py http://localhost:8080/ username hp.obo benchmark_data.tsv phenotips_results.tsv
    ```

#### Vibe

1. Follow steps for [running vibe][vibe_preperations].

2. Generate the bash scripts:
    ```
    python3 VibeBenchmarkParallelBashScriptsGenerator.py hp.obo benchmark_data.tsv ./ <MAX RUNS PER BASH SCRIPT>
    ```

3. Move the `vibe-with-dependencies.jar`, `TDB` and the bash scripts generated in the previous step
    (`vibe_benchmark_x.sh` where "x" indicates a number) to a folder where the benchmark will be done.

4. Rename/copy the `TDB` so that there are an equal amount of TDBs as there are benchmark bash scripts:
    ```
    # Rename the first TDB to be coherent with parallel bash script requirements.
    mv TDB/ TDB0/
    # Create copies of TDB equal to the number of created bash scripts (might require less/more copies than below).
    cp TDB0/ TDB1/
    cp TDB0/ TDB2/
    ```

5. Make final preparations for parallel benchmarking:
    ```
    mkdir results
    mkdir out
    mkdir err
    ```

6. Run benchmark:
    ```
    # Run all generated bash scripts (there might be less/more scripts than in this example).
    sh vibe_benchmark_0.sh
    sh vibe_benchmark_1.sh
    sh vibe_benchmark_2.sh
    ```

7. Process benchmark output:
    ```
    python3 VibeBenchmarkFileGenerator.py results/ vibe_results.tsv
    ```



[vibe]:https://github.com/molgenis/vibe
[vibe_preperations]:https://github.com/molgenis/vibe/#preparations
[hgnc_complete]:http://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
[query_phenomizer]:https://github.com/svandenhoek/query_phenomizer
[phenotips_download]:https://phenotips.org/Download

[hpo_obo_current]:http://purl.obolibrary.org/obo/hp.obo
[hpo_obo]:https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/2f6309173883d5d342849388c74bd986a2c0092c/hp.obo