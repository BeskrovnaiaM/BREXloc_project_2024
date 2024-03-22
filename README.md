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

### Repo content description

**scripts**

- `NCBI_data` - download small data *via* `ncbi-datasets`
- `NCBI_large_data` - download large data (> 1000 genomes or > 15 Gb) *via* `ncbi-datasets`
- `padloc_from_list` - runnig `padloc` on `ncbi-dataset` downloaded data for list of accessions IDs
- `padloc_from_folder` - runnig `padloc` on `ncbi-dataset` downloaded data for all subdirectories
- `extract_hmm_profiles` - extract hmm-profiles for Brx/Plg proteins from padloc hmm-database file
- `multibrex_search.py` - search genomes with more than one copy of BREX-system in `padloc` results
- `brex-extraction-padloc.ipynb` - search coordinates of BREX-system loci in `padloc` results, adds 10 kb flanks and creates corresponding proteins-fasta, gff-files and modified gff with proteins localisation label relative BREX-loci: upstream, inner or downstream.
- `non_brex_regs_search` - search non-brex systems on `padloc` results with additional gff-file with localisation label, which allows localization to be taken into account

**data**

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

**Links**
- [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/) CLI tool v16.4.5
- [Padloc](https://github.com/padloc/padloc) v2.0.0
