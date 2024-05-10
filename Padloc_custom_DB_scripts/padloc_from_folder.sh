#!/bin/bash
#v1.2
# Выполняет поиск по всем подкаталогам
# Для нуклеотидных последовательностей без аннотации
# Специализирован под результаты загрузки с помощью ncbi-dataset
# Аргументы:
###1 - путь, по которому лежат результаты ncbi-dataset
###2 - название папки для результатов
###3 - количество потоков

# Пример ввода команды:
# bash padloc_from_list.sh Downloads Results


mkdir $2

for VAR in $(ls $1/ncbi_dataset/data/ | grep '^GCF')
do
padloc --fna $1/ncbi_dataset/data/${VAR}/${VAR}*.fna --outdir $2 --cpu $3
done

