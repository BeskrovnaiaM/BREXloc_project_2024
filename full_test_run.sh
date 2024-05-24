#!/usr/bin/bash
# Large files (Pfam-A database) need to be downloaded
# Time consuming execution

wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_35.0.tar.gz

mkdir pfam

tar -xzf pfamA_35.0.tar.gz -C pfam

echo -e 'Downloading is complete\n\n'

tar -xzf test_data/test_padloc_data.tar.gz

python3 src/brex-extraction-padloc.py test_padloc_data BREX_loci BREX_loci_proteins.fasta BREX_loci_coordinates.txt

bash src/mmseqs_clusters.sh BREX_loci_proteins.fasta

bash src/hhsuite_annotation.sh unpacked_msa pfam/pfam

bash src/padloc_annotation.sh BREX_loci p_annotations

python3 src/annotation_table.py clustering_results/clusters_table.tsv \
p_annotations/ \
BREX_loci/ \
clustering_results/clusters_by_size.tsv \
hh_annotations/ \
test

