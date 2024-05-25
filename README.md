> **Bioinformatics Institute Course Project 2024**
# Search of novel BREX-associated immune systems with complementary mode of action, based on comparative genomics approach
### *Margarita Beskrovnaia, Kirill Petrikov*
### Supervised by Oksana Kotovskaya

&nbsp;  
**Aim**: To characterize the diversity of gene clusters within BREX loci and implement a pipeline to search this kind of loci in metagenomics data.

**Objectives:**
- To implement a workflow for extraction (M) and characterization of variable gene clusters within BREX loci
- To validate the workflow on the open databases
- To arrange the implemented workflow in the form of a Docker container to analyze newly sequenced data

### Content

- [Docker-meta for brexctractor](#Docker-meta)

- [Main workflow](#Main-workflow)

- [Test run](#Test-run)

- [Padloc_custom_DB folders](#Padloc_custom_DB-folders)

- [Links](#Links)

## `Docker-meta`

Contains the files needed to build a docker image that allows extraction of brex regions from a metagenome database.

**Usage**

Download files

```bash
mkdir brextracror
cd brextracror
wget https://raw.githubusercontent.com/BeskrovnaiaM/BREXloc_project_2024/progress/Docker-meta/Dockerfile \
https://raw.githubusercontent.com/BeskrovnaiaM/BREXloc_project_2024/progress/Docker-meta/brextractor.py \
https://raw.githubusercontent.com/BeskrovnaiaM/BREXloc_project_2024/progress/Docker-meta/requirements.txt
```

Build a docker image

```bash
docker build -t brextractor:1 .
```

Run container, mount the directory with metagenomic data to `/home/brextraction`

```bash
docker run -it --rm --name brextractor --rm -v [PATH_TO_DATA]:/home/brextraction brextractor:1
```

In container create directory for output files and run the `brextractor.py`

```bash
mkdir brextraction/results
python3 brextractor.py brextraction/[DIR_GBK] brextraction/results
```

`brextractor.py` takes two positional arguments:
- path to the directory with metagenomic data in gbk-file format
- path to the directory for output files

The outputs are:
- `output.fasta`: all protein sequences from extracted regions
- `regions_gff`: directory with custom gff-file, one for each region

**Important**
- The directory for output files must first be created
- Each time the script is launched, it adds new entries to the file and does not overwrite it

The script processes files in a single directory at once. If you need to process e.g. all directories, you can do this using a simple loop

```bash
for FILE in $(ls brexctraction); do python3 brextractor.py brextraction/${FILE} brextraction/results/; done
```

## Main workflow
The scripts for the main workflow are in folder `src`

## `brex-extraction-padloc.py`
Searches for genomic regions containing BREX defense systems. Selects these regions with specified up- and downstrean flanks (10 Kb each). Adds localisation mark to each protein: `brex`, `inner`, `upstream` or `downstream`.

**Requires**

PADLOC search results for nucleotide sequences

**Input positional arguments**
- path to directory with PADLOC results
- path to output directory for BREX-loci files
- output `fasta`-file name for all proteins
- output `txt`-file name for list of BREX-loci coordinates

**Output**
- `fasta`-file with all proteins from the BREX-loci found in the given contigs
- `txt`-file with BREX-loci coordinates
- for each contig creates regions files at the specified output directory:
  - `[PREFIX]_reg.faa` - `fasta`-file with proteins from given BREX-loci
  - `[PREFIX]_reg_positions.gff` - `gff`-file with added localisation marks
  - `[PREFIX]_reg.gff` - `gff`-file for given BREX-loci

The `[PREFIX]` is composed as follows:

`[Genome NCBI accession ID]_[BREX system name]_[PADLOC system number]_[Nucleotide NCBI ID]`
> Such scheme allows you to take into account the presence of systems of the same type in one nucleotide

### `mmseqs_clusters.sh`

Performs clustering by `MMseqs2` tool with the following parameters: `--min-seq-id 0.30 -c 0.8 --cov-mode 1`

The clusters names are composed as name of representative protein with adding `CLUS_` prefix

**Requires**
- `fasta`-file with all proteins from the BREX-loci (`brex-extraction-padloc.py` output)

**Input positional arguments**
- path to `fasta`-file

**Main output**
Create `clustering_results` directory with, among other:
- 'clusters_table.tsv' - table `Cluster-Protein`
- 'clusters_by_size.tsv' - table `Cluster-Cluster size`

Create `unpacked_msa/` directory with multiple sequence alignments for each cluster

> All logs are written to a `mmseqs.log` 

### `hhsuite_annotation.sh`

Performs annotatin by `HH-Suit` tool using `hhblits` command with one iteration (parameter `-n 1`)

**Requires**
- multiple sequence alignments for clusters (`mmseqs_clusters.sh` output)
- database for annotation

**Important**
The database must be preloaded

**Input positional arguments**
- path to directory with multiple sequence alignments
- path to database

**Output**
Create `hh_annotations` directory with annotations 'hhr'-files for each cluster

> All logs are written to a `hh_annotation.log`

### `padloc_annotation.sh`

Performs annotation by `PADLOC` tool

**Requires**
- directory with BREX-loci regions files (`brex-extraction-padloc.py` output)

**Input positional arguments**
- path to directory with BREX-loci regions files
- path to output directory

**Output**
- Standart PADLOC output files at the specified directory

> All logs are written to a `padloc_annotation.log`

### `annotation_table.py`

Creates summary tables with a list of clusters, clusters size and annotations, proteins and them localization marks.

**Requires**
- table `Cluster-Protein` (`mmseqs_clusters.sh` output)
- `PADLOC` annotations (`padloc_annotation.sh` output)
- directory with BREX-loci regions files (`brex-extraction-padloc.py` output)
- table `Cluster-Cluster size` (`mmseqs_clusters.sh` output)
- `HH-Suit` annotations (`hhsuite_annotation.sh` output)

**Input positional arguments**
- path to all requirements in the specified order
- prefix for output tables

**Output**
- `[PREFIX]_final_clusters_table.tsv` - table with selected single annotation for each cluster 
- `[PREFIX]_all_ann_clusters_table.tsv` - table with both `PADLOC` and `HH-Suit` annotations 

### `loci_counter.py`

Counts the number of regions with a unique composition in accordance with the given annotations. For a given cluster size threshold will count corresponding genes as identical thus loci with them as non unique. Such proteins will be renamed as `singleton`.

**Requires**
- table with selected single annotation for each cluster (`annotation_table.py` output)

**Input positional arguments**
- path to table
- cluster size to set threshold
- output `txt`-file name

**Main output**
- `txt`-file with list unique loci (by compositions as protein annotations according to gene location in the locus) and abundance

## `Test run`

The test scripts uses as input `test_padloc_data` - data obtained from the annotation of genomes by `PADLOC` with a custom database (only BREX-proteins hmm-profiles)

You can test run the workflow in two ways:
- `short_test_run.sh` uses ready-made clusterin results by `MMseqs2` and annotations results by `PADLOC` and `HH-Suite` (`clustering_results_ready`, `p_annotations_ready`, `hh_annotations_ready`).
- `full_test_run.sh` performs the full workflow, but is time-consuming as it requires downloading the `Pham-A` database

If all steps are successfully completed, two files will be created in the root directory, among other:
- `test_final_clusters_table.tsv`
- `test_loci_results.txt`

Examples of these expected results can also be found in `test_data`

> Note that clustering results may differ between different runs on the same data, so small amounts are acceptable.

**Instruction**

- Clone repo
```bash
git@github.com:BeskrovnaiaM/BREXloc_project_2024.git
```

- Create new environment `BREX_loci`
```bash
mamba env create -f environment.yml
```

- Run `short_test_run.sh` or `full_test_run.sh`
```bash
bash short_test_run.sh
bash full_test_run.sh
```

## `Padloc_custom_DB` folders

Contains materials for working with a test custom database

**additional_src**
- `NCBI_large_data.sh` - download large data (> 1000 genomes or > 15 Gb) *via* `ncbi-datasets`
- `extract_hmm_profiles` - extract hmm-profiles for Brx/Plg proteins from `PADLOC` hmm-database file
- `padloc_custom.sh` - runnig `PADLOC` on data downloaded by`NCBI_large_data.sh`. Uses custom hmm-database (only Brx/Plg proteins). **Requires** `extract_hmm_profiles.sh`
- `multibrex_search.py` - search genomes with more than one copy of BREX-system in `PADLOC` results
- `non_brex_regs_search` - search non-brex systems on `PADLOC` results with additional `gff`-file with localisation label, which allows localization to be taken into account

**Padloc_trial_DB_data**

- `2024_01_df-RefSeq` - data of bacteria containing BREX-systems.
  >Source: [DefenseFinder RefSeq DB](https://defensefinder.mdmlab.fr/wiki/refseq/) (date of access 31.01.2024)
- `2024_02_multi_brex_df-RefSeq` - entries from `2024_01_df-RefSeq` corresponding bacteria containing multiple BREX-systems. Column 'abundance' was added.
- `2024_02_genomes_accessios_dedup` - RefSeq accession IDs, extracted from `2024_01_df-RefSeq` after duplicate removal
- `2024_02_brex_hmm_list` - hmm-profiles for BREX-proteins extrated from `PADLOC` hmm-file *via* `extract_hmm_profiles`
- `2024_02_multi_brex_Padloc` - list of accessions with more than one copy of BREX-system *via* `multibrex_search`
- `2024_03_19_brex_regs_coord.bed` - BREX-regions coordinates

[**Internal data**](https://figshare.com/s/643c4203c7d2769bb938)

- `2024-02-26_padloc_hmm_custom` - `PADLOC` results for `2024_02_genomes_accessios_dedup` accessions (exept two entries: GCF_016887565.1, GCF_005931095.1)
- `2024_03_20_brex_extr_regs` - BREX-system loci with *c.* 10 kb flanks generated by `brex-extraction-padloc.py`. Contains corresponding `fasta-proteins` and `gff`-files, and modified `gff` with localisation column.
- `2024_03_20_prot_brex_regs` - all proteins from BREX-loci

## Links
- [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/) CLI tool v16.4.5
- [DefenseFinder RefSeq DB](https://defensefinder.mdmlab.fr/wiki/refseq) Database relevance: 05/2022
- [Padloc](https://github.com/padloc/padloc) v2.0.0
- [MMseqs2](https://github.com/soedinglab/MMseqs2) Release 15-6f452
- [HH-Suite](https://github.com/soedinglab/hh-suite) v3.3.0
