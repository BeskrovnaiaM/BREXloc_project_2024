"""v1.1
    Считает общее число не-brex-систем в csv-файлах Padloc,
    учитывая локализацию ЗС относительно brex-системы
    
    Первый аргумент - путь до директории с результатами Падлока
    Второй аргумент - путь до директории с gff-файлами разметки
    
    Записывает три tsv-файла:
    'Total_non_brex_systems.tsv' - все не-brex-системы
    'Total_inner_systems.tsv' - только с внутренней локализацией
    'Total_other_systems.tsv' - остальные
"""

import os
import sys
import pandas as pd


def write_results(output_file_name: str, resuls_to_write: dict) -> None:
    """
    Принимает словарь с подсчитанными ЗС и записывает в файл
    """
    with open(output_file_name, mode='w') as file:
        file.write('System_type' + '\t' + 'Number' + '\n')
        for key, value in resuls_to_write.items():
            new_line = key + '\t' + str(value)  +'\n'
            file.write(new_line)


# Первый аргумент - путь до директории с результатами Падлока
input_padloc_dir = sys.argv[1]

# Второй аргумент - путь до директории с gff-файлами, содержащими метку локализации 'position'
input_pos_dir = sys.argv[2]

# Словари для итогового подсчёта ЗС
total_non_brex_def_systems = {}  # суммарных не-brex-систем
total_inner_def_systems = {}  # только со внутренней локализацией
total_other_def_systems = {}  # с прочей локализацией

loci_with_inner_sys = 0  # счётчик регионов, где в принципе встречаются ЗС с внутренней локализацией

# Заголовок для gff-датафрейма
gff_header = ['seqid', '_2', '_3', '_4', '_5', '_6', '_7', '_8', 'Comment', 'Loc']

for file_name in os.listdir(input_padloc_dir):

    curr_non_brex_def_systems = {}
    curr_inner_def_systems = {}
    curr_other_def_systems = {}

    if file_name.endswith('csv'):
        # Таблица с результатами Padloc
        df = pd.read_csv(os.path.join(input_padloc_dir, file_name)).iloc[:, :4]

        # Таблица с меткой локализации белков
        file_name_pos = file_name.split('padloc')[0] + 'positions.gff'
        df_pos = pd.read_csv(os.path.join(input_pos_dir, file_name_pos),
                             sep='\t', names=gff_header).iloc[:, [0, 8, 9]]

        # создание уникального ID позволяет учесть совпадающие названия типов систем
        df['Name'] = df['system'] + '%' + df['system.number'].astype('str')

        # Столбик 'target.name' получается как и аналогичный в 'df'
        # Это позволит провести сопоставление локализации белков с принадлежностью к конкретной ЗС
        df_pos['target.name'] = (df_pos['seqid'] +
                                 '_' +
                                 [i for i in map(lambda x: x.split(';')[0].split('_')[1], df_pos['Comment'])])

        # Слияние получается по колонке'target.name', для не-брекс-систем, и из 'df_pos' добавляется локализация
        df_non_brex = df[df['system'].str[:4] != 'brex'].merge(df_pos.iloc[:, [2, 3]], how='left')

        # Группировка по уникальному для каждой системы 'Name', и создаётся соответствующая ячейка с локализацией
        uniq_sys = df_non_brex.loc[:, ('Name', 'Loc')].groupby('Name').agg(lambda x: set(x))

        # Список со всеми не-брекс-системами
        curr_non_brex_def_systems = [i for i in map(lambda x: x.split('%')[0], uniq_sys.index.unique())]
        
        # Только с внутренней локализацией
        curr_inner_def_systems = [i for i in
                                  map(lambda x: x.split('%')[0], uniq_sys[uniq_sys['Loc'] == {'inner'}].index.unique())]
        
        # Остальные
        curr_other_def_systems = [i for i in
                                  map(lambda x: x.split('%')[0], uniq_sys[uniq_sys['Loc'] != {'inner'}].index.unique())]
                
        if len(curr_inner_def_systems) > 0:
            loci_with_inner_sys += 1
                    
        # Обновляем подсчёт не-брекс-систем
        for def_sys in curr_non_brex_def_systems:
            if def_sys in total_non_brex_def_systems:
                total_non_brex_def_systems[def_sys] += 1
            else:
                total_non_brex_def_systems[def_sys] = 1

        for def_sys in curr_inner_def_systems:
            if def_sys in total_inner_def_systems:
                total_inner_def_systems[def_sys] += 1
            else:
                total_inner_def_systems[def_sys] = 1

        for def_sys in curr_other_def_systems:
            if def_sys in total_other_def_systems:
                total_other_def_systems[def_sys] += 1
            else:
                total_other_def_systems[def_sys] = 1

#  Сохраняем результаты
output_files_names = ('Total_non_brex_systems.tsv', 'Total_inner_systems.tsv', 'Total_other_systems.tsv')
def_systems = (total_non_brex_def_systems, total_inner_def_systems, total_other_def_systems)

for file_name, sys_dict in zip(output_files_names, def_systems):
    write_results(file_name, sys_dict)

