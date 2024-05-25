#!/bin/bash

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

