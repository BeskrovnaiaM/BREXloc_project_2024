"""v2.1
   Проверяет уникальность brex-систем по csv-файлам Padloc.
   Сохраняет список всех мульти-brex файлов в таблицу 'multi_brex_files.tsv'
   
   Может создать копию / перезаписать csv-файлы результатов Padloc,
   дополнив столбцом с уникальным ID вида: system_system.number_seqid
   
   Аргументы
   1 - список названий csv-файлов Padloc
       Предполагается, что все файлы лежат в одной директориию
       Имена файлов записываются в итоговую таблицу со списком мульти-brex файлов
   2 - путь до папки с этими файлами
   3 необязательный - указать имя папки, если нужно сохранить csv-файлы Padloc с дополненным столбцом,
                      если указать REWRITE - добавит столбец, перезаписав исходные csv-файлы
"""

import os
import sys
import pandas as pd


def myf(file_path, file_name, write_csv=False):
    """
    Из столбца 'system' отбирает только 'brex'.
    Считает, сколько уникальных 'system.number'.
    Если больше 1 - возвращает имя файла, количество brex-систем,
                    количество контиг для brex,
                    и список brex-систем с именами вида:
                    system_system.number_seqid, например:
                    brex_type_I_1_NC_008346.1
    
    Сохраняет csv-файл с добавляя в конец столбец 'Name' с именами,
    если указан путь write_csv
    """

    df = pd.read_csv(file_path)

    df['Name'] = df['system'] + '_' + df['system.number'].astype('str') + '_' + df['seqid']

    if write_csv:
        path_to_save = os.path.join(write_csv, file_name)
        df.to_csv(path_to_save, index=False)

    brex_contigs_number = len(df[df['system'].str.startswith('brex')]['seqid'].unique())

    uniq_brex_number = len(df[df['Name'].str.startswith('brex')]['system.number'].unique())

    uniq_brex_names = df[df['Name'].str.startswith('brex')]['Name'].unique()

    if uniq_brex_number > 1:
        return file_name, str(uniq_brex_number), str(brex_contigs_number), uniq_brex_names.tolist()


# Первый аргумент - список csv-файлов Падлока
input_files_list = sys.argv[1]

# Второй аргумент - путь до csv-файлов
input_padloc_dir = sys.argv[2]

# Третий аргумент необязательный.
# Если указан, то это будет имя директории для сохранения дополненных csv-файлов
# Если REWRITE - перезапишет исходные csv-файлы
if len(sys.argv) == 4:
    if sys.argv[3] == 'REWRITE':
        output_dir = input_padloc_dir
    else:
        output_dir = sys.argv[3]
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
else:
    output_dir = False

# Прочитать список с названиями файлов для сканирования
files_list = []

with open(input_files_list) as file:
    for line in file:
        if not line.strip() == '':  # Проверка на пустые строки
            files_list.append(line.strip())

# Найти мульти-брекс системы
list_multi_brex_files = []

for name in files_list:
    res = myf(os.path.join(input_padloc_dir, name), name, output_dir)
    if res:
        list_multi_brex_files.append(res)

# Запись результатов
header = 'Accession' + '\t' + 'Total_brex' + '\t' + 'Contigs_brex' + '\t' + 'Brex_list' + '\n'
with open(os.path.join('multi_brex_files.tsv'), mode='w') as file:
    file.write(header)
    for record in list_multi_brex_files:
        f_name, b_num, c_num, b_names = '.'.join(record[0].split('.')[:-1]), record[1], record[2], ','.join(record[3])
        new_rec = f_name + '\t' + b_num + '\t' + c_num + '\t' + b_names + '\n'
        file.write(new_rec)

