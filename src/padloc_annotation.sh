#!/usr/bin/bash
# Args:
###1 path to folder with input data
###2 path to output folder


mkdir $2

echo 'Start Padloc'

for VAR in $(ls $1 | grep 'faa')
do
 NAME=${VAR::-4}
 echo "Processing ${NAME}"
 
 padloc --faa $1/${NAME}.faa --gff $1/${NAME}.gff --outdir $2 --fix-prodigal >>padloc_annotation.log 2>&1
done

echo 'Finish Padloc'

