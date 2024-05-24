#!/usr/bin/bash
# Arg:
###1 path to msa-files
###2 path to pfam hh-suite database

mkdir hh_annotations

echo 'Start HH-Suite'

for FILE in $(ls $1)
do
 echo "Processing ${FILE}"
 hhblits -n 1 -i $1/${FILE} -d $2 -o hh_annotations/${FILE}.hhr >>hh_annotation.log 2>&1
done

echo 'Finish HH-Suite'

