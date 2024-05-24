#!/usr/bin/bash
# Arg:
###1 path to fasta-file

LINE='-------------------\n'


echo -e "${LINE}mmseqs2 Protein clustering\n${LINE}"
echo -e "${LINE}mmseqs2 Protein clustering\n${LINE}" >>mmseqs.log 2>&1

mmseqs createdb $1 mmseqs_db >>mmseqs.log 2>&1
mmseqs cluster mmseqs_db mmseqs_clusters_db mmseqs_tmp --min-seq-id 0.30 -c 0.8 --cov-mode 1 >>mmseqs.log 2>&1


echo -e "${LINE}mmseqs2 Create msa\n${LINE}"
echo -e "${LINE}mmseqs2 Create msa\n${LINE}" >>mmseqs.log 2>&1

mmseqs result2msa mmseqs_db mmseqs_db mmseqs_clusters_db msa_clusters --msa-format-mode 5 >>mmseqs.log 2>&1
mkdir unpacked_msa
mmseqs unpackdb msa_clusters unpacked_msa/ >>mmseqs.log 2>&1


echo -e "${LINE}Rename files\n${LINE}"

for FILE in $(ls unpacked_msa/)
do
 NAME="CLUS_$(head -1 unpacked_msa/${FILE} | sed 's/^.//')"
 mv unpacked_msa/${FILE} unpacked_msa/${NAME}
done
echo -e "${LINE}Renaming complete\n${LINE}"


echo -e "${LINE}mmseqs2 Create clustering results\n${LINE}"
echo -e "${LINE}mmseqs2 Create clustering results\n${LINE}" >>mmseqs.log 2>&1

mkdir clustering_results/

mmseqs createtsv mmseqs_db mmseqs_db mmseqs_clusters_db clustering_results/clusters.tsv >>mmseqs.log 2>&1

### Создаю фасту
mmseqs createseqfiledb mmseqs_db mmseqs_clusters_db fasta_db >>mmseqs.log 2>&1

mmseqs result2flat mmseqs_db mmseqs_db fasta_db clustering_results/proteins_by_clusters.fasta >>mmseqs.log 2>&1

### Создать список представителей кластера (центроидов). Сначала база данных
mmseqs createsubdb mmseqs_clusters_db mmseqs_db representatives_db >>mmseqs.log 2>&1

mmseqs convert2fasta representatives_db clustering_results/representatives.fasta >>mmseqs.log 2>&1

echo -e "${LINE}mmseqs2 Create clustering summary\n${LINE}"

awk 'BEGIN{FS="\t";OFS="\t"} {print "CLUS_"$1,$2}' clustering_results/clusters.tsv >clustering_results/clusters_table.tsv
cat clustering_results/clusters.tsv | cut -f 1 | sort | uniq -c | awk 'BEGIN{FS=" ";OFS="\t"} {print "CLUS_"$2,$1}' >clustering_results/clusters_by_size.tsv


rm -rf mmseqs_tmp

