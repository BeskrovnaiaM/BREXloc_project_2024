#!/bin/bash
#v1.1
# Для загрузки из Генбанка больших объёмов данных: > 1000 геномов или > 15 Гб (dehydrated/rehydrate способ)
# Переименовывает файлы 'genomic.gff' и 'protein.faa', добавляя к названию соответствующий ID сборки
# Два аргумента:
### Первый - название файла со списком ID сборок
### Второй - название папки, куда сохранятся результаты
# Создаёт файл 'Processed_IDs_list.txt' с ID обработанных сборок

# Пример ввода команды:
# bash NCBI_large_data.sh accessions_list.txt Genomes


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

# Логирование списка загруженных геномов
echo ${NAME} >>Processed_IDs_list.txt
done

echo "---Files renamed. 'Processed_IDs_list.txt' created---"

# Удаление служебных файлов
rm $2.zip
rm $2/ncbi_dataset/fetch.txt
rm $2/README.md
