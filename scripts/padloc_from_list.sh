#!/bin/bash
#v1.2
# Выполняет поиск по ID из файла
# Для нуклеотидных последовательностей без аннотации
# Специализирован под результаты загрузки с помощью ncbi-dataset
# Аргументы:
###1 - путь, по которому лежит список ID
###2 - путь, по которому лежат результаты ncbi-dataset
###3 - название папки для результатов

# Пример ввода команды:
# bash padloc_from_list.sh list.txt Downloads Results


mkdir $3

for VAR in $(cat $1)
do
padloc --fna $2/ncbi_dataset/data/${VAR}/${VAR}*.fna --outdir $3
done

