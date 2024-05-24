#!/usr/bin/bash

tar -xzf test_data/test_padloc_data.tar.gz
tar -xzf test_data/p_annotations_ready.tar.gz
tar -xzf test_data/hh_annotations_ready.tar.gz
tar -xzf test_data/clustering_results_ready.tar.gz

python3 src/brex-extraction-padloc.py test_padloc_data BREX_loci BREX_loci_proteins.fasta BREX_loci_coordinates.txt

python3 src/annotation_table.py clustering_results_ready/clusters_table.tsv \
p_annotations_ready/ \
BREX_loci/ \
clustering_results_ready/clusters_by_size.tsv \
hh_annotations_ready/ \
test

python3 src/loci_counter.py test_final_clusters_table.tsv 0 test_loci_results.txt

