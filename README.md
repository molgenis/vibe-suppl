# vibe-suppl
This repo contains supplemental files regarding the Java application found [here][vibe]. Note that these are in no way needed to use the vibe tool, but were used to generate additional information (such as benchmarking). They were created with the assumption that they are used exactly in the way they are meant to be used, so while certain checks/validations might be present, using these scripts in the wrong way might result in weird behavior.

## Paper

Please refer to the `README.md` at https://zenodo.org/record/3662470 for the exact commits used for the benchmarking. There, all required files for [PaperPlots.R](benchmarking_results_processing/PaperPlots.R) can be found as well.

## Benchmarking

### Data

There are several files used among these scripts. These include:
* [benchmark_data.tsv](https://zenodo.org/record/3662470/files/benchmark_data-hgnc_symbol.tsv)
    * A dataset with the first column being an ID and the fourth column 1 or more phenotypes separated
    by a comma (the phenotype names should exist within the [Human Phenotype Ontology][hpo_obo]) .
* [hp.obo][hpo_obo]
    * The Human Phenotype Ontology used for combining/converting phenotype names with their HPO ID. Note that the `benchmark_data.tsv` was made compatible for release 2018-03-08 specifically.
* [hgnc_complete_set.txt][hgnc_complete]
    * The HUGO Gene Nomenclature Committee file containing information about genes (primarily used to generate a list containing all genes).
* [benchmark_file_conversion_data.tsv](https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_prev_sym&col=md_eg_id&col=gd_pub_eg_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_hgnc_id&format=text&submit=submit)
  * A file generated through [genenames.org](https://www.genenames.org/) that contains HGNC gene symbols with their previous symbols and their NCBI gene IDs.

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

3. Convert the HGNC gene symbols to NCBI gene IDS:

    ```
    python3 BenchmarkFileGeneSymbolToIdConverter.py amelie_results.tsv benchmark_file_conversion_data.tsv 1> amelie.log 2> amelie.err
    ```

#### Exomiser

**IMPORTANT:** A custom `.jar` file supplied by the Exomiser team was supplied to run this benchmark without requiring a `.vcf` file. Exomiser has not yet made a public release of this yet. This custom `.jar` however is based on the exomiser-rest-prioritiser module of the Exomiser open-source code (release 12.1.0).

##### hiPHIVE

1. Run benchmark:

   ```
   python3 ExomiserBenchmarkRunner.py hp.obo benchmark_data.tsv hiphive hiphive_output/
   ```

2. Process benchmark output:

   ```
   python3 ExomiserBenchmarkFileGenerator.py hiphive_output/ hiphive_results.tsv
   ```

3. Convert the HGNC gene symbols to NCBI gene IDS:

   ```
   python3 BenchmarkFileGeneSymbolToIdConverter.py hiphive_results.tsv benchmark_file_conversion_data.tsv 1> hiphive.log 2> hiphive.err
   ```

##### PhenIX

1. Run benchmark:

   ```
   python3 ExomiserBenchmarkRunner.py hp.obo benchmark_data.tsv phenix phenix_output/
   ```

2. Process benchmark output:

   ```
   python3 ExomiserBenchmarkFileGenerator.py phenix_output/ phenix_results.tsv
   ```

3. Convert the HGNC gene symbols to NCBI gene IDS:

   ```
   python3 BenchmarkFileGeneSymbolToIdConverter.py phenix_results.tsv benchmark_file_conversion_data.tsv 1> phenix.log 2> phenix.err
   ```

#### GADO

We used the stand-alone commandline version GADO (v 1.0.1), available at: https://github.com/molgenis/systemsgenetics/wiki/GADO-Command-line. We accepted all automatically suggested alternative HPO terms in cases that the supplied HPO term could not be used. We have used the prediction matrix `hpo_predictions_sigOnly_spiked_01_02_2018`. The output was also converted to NCBI gene IDs through the following:

```
python3 BenchmarkFileGeneSymbolToIdConverter.py gado_results.tsv benchmark_file_conversion_data.tsv 1> gado.log 2> gado.err
```

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
    
4. Convert the HGNC gene symbols to NCBI gene IDS:

    ```
    python3 BenchmarkFileGeneSymbolToIdConverter.py phenomizer_results.tsv benchmark_file_conversion_data.tsv 1> phenomizer.log 2> phenimozer.err
    ```

#### Phenotips

**IMPORTANT**: As of January 2020, Phenotips does not offer a stand-alone downloadable solution anymore and requires a paid cloud subscription to be used ([source](https://phenotips.com/blog/new-year-new-website.html)). While the [GitHub repo](https://github.com/phenotips/phenotips) is currently still online, it seems uncertain whether it will still be updated and the easy-to-use `.dmg` as offered on the old website is not available anymore. Therefore, this benchmark is deemed obsolete.

#### PubCaseFinder

1. Run benchmark:

   ```
   python3 PubCaseFinderBenchmarkRunner.py hp.obo benchmark_data.tsv pubcasefinder_output/
   ```

2. Process benchmark output:

   ```
   python3 PubCaseFinderBenchmarkFileGenerator.py pubcasefinder_output/ pubcasefinder_results.tsv
   ```

3. Convert the HGNC gene symbols to NCBI gene IDS:

   ```
   python3 BenchmarkFileGeneSymbolToIdConverter.py amelie_results.tsv benchmark_file_conversion_data.tsv 1> amelie.log 2> amelie.err
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
    python3 VibeBenchmarkFileGenerator.py results/ vibe_results.tsv none
    ```

[vibe]:https://github.com/molgenis/vibe
[vibe_preperations]:https://github.com/molgenis/vibe/#quickstart
[hgnc_complete]:http://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
[query_phenomizer]:https://github.com/svandenhoek/query_phenomizer
[phenotips_download]:https://phenotips.org/Download

[hpo_obo_current]:http://purl.obolibrary.org/obo/hp.obo
[hpo_obo]:https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/2f6309173883d5d342849388c74bd986a2c0092c/hp.obo

