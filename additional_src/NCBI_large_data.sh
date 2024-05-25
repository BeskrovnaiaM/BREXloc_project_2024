#!/bin/bash
# To download large amounts of data from the GenBank: > 1000 genomes or > 15 GB (dehydrated/rehydrated method)
# Renames the files 'genomic.gff' and 'protein.faa', adding the appropriate assembly ID to the name
# Args:
###1 file with a list of assembly IDs to download (Each ID must be on a separate line)
###2 path to directory where the results will be saved

datasets download genome accession --inputfile $1 --dehydrated --include genome,gff3,protein --filename $2.zip

echo "---Dehydrated data download complete---"

unzip -q $2.zip -d $2

echo "---Files unpacked---"

datasets rehydrate --directory $2

echo "---Rehydration complete---"

for NAME in $(ls $2/ncbi_dataset/data/ | grep '^GCF')
do
mv $(ls $2/ncbi_dataset/data/${NAME}/*.gff) $2/ncbi_dataset/data/${NAME}/${NAME}_genomic.gff
mv $(ls $2/ncbi_dataset/data/${NAME}/*.faa) $2/ncbi_dataset/data/${NAME}/${NAME}_protein.faa
done

echo "---Files renamed. 'Processed_IDs_list.txt' created---"

rm $2.zip
rm $2/ncbi_dataset/fetch.txt
rm $2/README.md

