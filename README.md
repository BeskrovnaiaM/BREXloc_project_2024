> **Bioinformatics Institute Course Project 2024**
# Search of novel BREX-associated immune systems with complementary mode of action, based on comparative genomics approach
### *Margarita Beskrovnaia, Kirill Petrikov*
### Supervised by Oksana Kotovskaya
> Branch `progress` for the working version of the project

&nbsp;  
**Aim**: To characterize the diversity of gene clusters within BREX loci and implement a pipeline to search this kind of loci in metagenomics data.

**Objectives:**
- To implement a workflow for extraction (M) and characterization of variable gene clusters within BREX loci
- To validate the workflow on the open databases
- To arrange the implemented workflow in the form of a Docker container to analyze newly sequenced data

**Content**

[Docker-meta for brexctractor](#Docker-meta-for-brexctractor)

[Padloc_custom_DB folders](#Padloc_custom_DB-folders)

[Links](#Links)

## `Docker-meta` for `brexctractor`

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


## `Padloc_custom_DB` folders

Contains materials for working with a test custom database

**Padloc_custom_DB_scrips**

- `NCBI_data` - download small data *via* `ncbi-datasets`
- `NCBI_large_data` - download large data (> 1000 genomes or > 15 Gb) *via* `ncbi-datasets`
- `padloc_from_list` - runnig `Padloc` on `ncbi-dataset` downloaded data for list of accessions IDs
- `padloc_from_folder` - runnig `Padloc` on `ncbi-dataset` downloaded data for all subdirectories
- `extract_hmm_profiles` - extract hmm-profiles for Brx/Plg proteins from padloc hmm-database file
- `multibrex_search.py` - search genomes with more than one copy of BREX-system in `Padloc` results
- `brex-extraction-padloc.ipynb` - search coordinates of BREX-system loci in `Padloc` results, adds 10 kb flanks and creates corresponding proteins-fasta, gff-files and modified gff with proteins localisation label relative BREX-loci: upstream, inner or downstream.
- `non_brex_regs_search` - search non-brex systems on `Padloc` results with additional gff-file with localisation label, which allows localization to be taken into account

**Padloc_custom_DB_data**

- `2024_01_df-RefSeq` - data of bacteria containing BREX-systems.
  >Source: [DefenseFinder RefSeq DB](https://defensefinder.mdmlab.fr/wiki/refseq/) (date of access 31.01.2024)
- `2024_02_multi_brex_df-RefSeq` - entries from `2024_01_df-RefSeq` corresponding bacteria containing multiple BREX-systems. Column 'abundance' was added.
- `2024_02_genomes_accessios_dedup` - RefSeq accession IDs, extracted from `2024_01_df-RefSeq` after duplicate removal
- `2024_02_brex_hmm_list` - hmm-profiles for BREX-proteins extrated from `Padloc` hmm-file *via* `extract_hmm_profiles`
- `2024_02_multi_brex_Padloc` - list of accessions with more than one copy of BREX-system *via* `multibrex_search`
- `2024_03_19_brex_regs_coord.bed` - BREX-regions coordinates

[**Internal data**](https://figshare.com/s/643c4203c7d2769bb938)

- `2024-02-26_padloc_hmm_custom` - `padlock` results for `2024_02_genomes_accessios_dedup` accessions (exept two entries: GCF_016887565.1, GCF_005931095.1)
- `2024_03_20_brex_extr_regs` - BREX-system loci with *c.* 10 kb flanks extracted *via* `brex-extraction-padloc.ipynb`. Contains corresponding proteins-fasta and gff-files, and modified gff with localisation label.
- `2024_03_20_prot_brex_regs` - all proteins from BREX-regions

## Links
- [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/) CLI tool v16.4.5
- [Padloc](https://github.com/padloc/padloc) v2.0.0
