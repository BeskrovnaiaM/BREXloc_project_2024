#!/bin/bash
#v1.2
# Search loc by padloc with proteins-fasta and gff-annotation
# Специализирован под результаты загрузки с помощью ncbi-dataset
# Принимает файлы:
### аминокислотные последовательности .faa; разметку .gff
# Файлы последовательностей должны находиться в /ncbi_dataset/data/, каждый в своей папке с тем же названием, что и сам файл
# Три аргумента:
### Первый - путь, по которому лежит папка ncbi-dataset/
### Второй - название папки для результатов
### Третий - число потоков

# Пример ввода команды:
# bash padloc_prot_script.sh Genomes Results 8


mkdir $2

for VAR in $(ls $1/ncbi_dataset/data/ | grep '^GCF')
do
padloc --faa $1/ncbi_dataset/data/${VAR}/${VAR}*.faa --gff $1/ncbi_dataset/data/${VAR}/${VAR}*.gff --outdir $2 --cpu $3
done

