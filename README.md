# BREXloc_project_2024

> Branch `progress` for the working version of the project

### Repo content description

**scripts**

- `NCBI_data` - download small data *via* `ncbi-datasets`
- `NCBI_large_data` - download large data (> 1000 genomes or > 15 Gb) *via* `ncbi-datasets`
- `padloc_from_list` - runnig `padloc` on `ncbi-dataset` downloaded data for list of accessions IDs
- `padloc_from_folder` - runnig `padloc` on `ncbi-dataset` downloaded data for all subdirectories
- `extract_hmm_profiles` - extract hmm-profiles for Brx/Plg proteins from padloc hmm-database file

**data**

- `2024_01_df-RefSeq` - list of BREX-system contained bacteria obtained from [DefenceFinder RefSeq DB](https://defensefinder.mdmlab.fr/wiki/refseq/) (date of access 31.01.2024)
- `2024_02_genomes_accessios_dedup` - RefSeq accession IDs, extracted from `2024_01_df-RefSeq` after duplicate removal
- `brex_hmm_list` - hmm-profiles for BREX-proteins extrated from `Padloc` hmm-file *via* `extract_hmm_profiles`.
