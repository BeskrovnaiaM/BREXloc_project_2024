#!/usr/bin/bash
# Args:
###1 path to padloc database
###2 path to genomes downloaded by NCBI Datasets

LINE='-------------------\n'

cp $1/padlocdb.hmm ./

python3 extract_padloc_brex_hmm.py $1

rm $1/padlocdb.hmm

cp brex_select.hmm $1/padlocdb.hmm

echo -e "${LINE}Hmm-profile updated to custom file\n${LINE}"

mkdir padloc_brex

for VAR in $(ls $2/ncbi_dataset/data/ | grep '^GCF')
do
 padloc --fna $2/ncbi_dataset/data/${VAR}/${VAR}*.fna --outdir padloc_brex/
done

echo -e "${LINE}Custom hmm-padloc complete\n${LINE}"

rm $1/padlocdb.hmm
cp padlocdb.hmm $1

echo -e "${LINE}Default hmm-profile restored\n${LINE}"

