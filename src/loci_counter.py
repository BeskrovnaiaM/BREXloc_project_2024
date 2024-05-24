import csv
import sys

from collections import Counter


def update_second_column(input_file, output_file):
    updated_lines = []
    with open(input_file, 'r') as f:
        for line in f:
            if line.strip():  # проверяем, что строка не пустая
                parts = line.split()
                if len(parts) >= 2:
                    second_column_value = parts[2]
                    if '_' in second_column_value:
                        prefix, suffix = second_column_value.rsplit('_', 1)
                        if len(suffix) < 5:
                            suffix = suffix.zfill(5)
                        updated_second_column_value = prefix + '_' + suffix
                        parts[2] = updated_second_column_value
                updated_line = ' '.join(parts) + '\n'
                updated_lines.append(updated_line)

    with open(output_file, 'w') as f:
        f.writelines(updated_lines)


def get_number_from_protein(protein_string):
    # Находим второе появление символа "_" и извлекаем число после него
    second_underscore_index = protein_string.find('_', protein_string.find('_') + 1)
    number_string = protein_string[second_underscore_index + 1:]
    # Преобразуем строку в число, учитывая ведущие нули
    return int(number_string.lstrip('0'))


def write_to_output_file(clust_region_list, output_file):
    with open(output_file, 'a') as f:
        f.write(','.join(clust_region_list) + '\n')


def main(input_file, output_file, small_cluster_threshold=0):
    # Сортируем файл по второй колонке
    with open(input_file, 'r') as f:
        lines = f.readlines()
        # Пропускаем первую строку (хэдер)
        header = lines[0]
        data = sorted(lines[1:], key=lambda x: x.split(' ')[2])

    clust_region = []
    prev_protein_number = None

    for line in data:
        # Разделяем строку по пробелу
        columns = line.strip().split(' ')
        # Проверяем, есть ли четыре столбца
        if len(columns) < 3:
            print("Ошибка: Неверное количество столбцов в строке:", line)
            continue
        protein = columns[2]
        size_clust = int(columns[1])
        clust_region_value = columns[0]
        if size_clust <= small_cluster_threshold:
            clust_region_value = 'x' + clust_region_value

        # Получаем число из второй колонки (protein)
        try:
            protein_number = get_number_from_protein(protein)
        except ValueError:
            print("Ошибка: Невозможно извлечь число из строки:", protein)
            continue

        if prev_protein_number is None:
            prev_protein_number = protein_number

        # Проверяем разницу между текущим и предыдущим числами
        if abs(protein_number - prev_protein_number) <= 1:
            # Добавляем значение в лист clust_region
            clust_region.append(clust_region_value + '_' + columns[-2])  # Взять значение из последней колонки
        else:
            # Разница больше 1, записываем текущий список в новый файл
            write_to_output_file(clust_region, output_file)
            # Очищаем список
            clust_region = []
            # Добавляем текущее значение в новый список
            clust_region.append(clust_region_value + '_' + columns[-2])  # Взять значение из последней колонки

        prev_protein_number = protein_number

    # Записываем оставшиеся данные в файл
    write_to_output_file(clust_region, output_file)


def reverse_proteins(line):
    proteins = line.strip().split(',')
    inner_index = None
    brex_index = None
    for i, protein in enumerate(proteins):
        if protein.endswith('_inner') or protein.endswith('_brex'):
            inner_index = i
            brex_index = i
            break
    if inner_index is not None:
        downstream_index = None
        for i, protein in enumerate(proteins):
            if protein.endswith('_downstream'):
                downstream_index = i
                break
        if downstream_index is not None and downstream_index < inner_index:
            reversed_proteins = proteins[::-1]
            return ','.join(reversed_proteins)
    return line.strip()


def filter_csv(input_file, output_file):
    filtered_rows = []

    with open(input_file, 'r', newline='') as f_in:
        reader = csv.reader(f_in)
        for row in reader:
            filtered_row = [item for item in row if item.endswith('_brex') or item.endswith('_inner')]
            filtered_rows.append(filtered_row)

    with open(output_file, 'w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerows(filtered_rows)


def process_csv(input_file, output_file):
    with open(input_file, 'r') as infile:
        with open(output_file, 'w') as outfile:
            for line in infile:
                elements = line.strip().split(',')
                modified_elements = [element if not element.startswith('x')
                                     else 'SINGLETONE' for element in elements]
                modified_line = ','.join(modified_elements) + '\n'
                outfile.write(modified_line)


def process_protein_names(names):
    # Удаляем приписки "_inner" и "_brex" и сортируем имена белков
    processed_names = [name.replace("_inner", "").replace("_brex", "") for name in names]
    return tuple(processed_names)


def count_unique_regions(file_path):
    regions_counter = Counter()

    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            # Преобразуем имена белков и считаем уникальные регионы
            region = process_protein_names(row)
            regions_counter[region] += 1

    # Сортируем словарь по значениям (количество встречаемости региона) в убывающем порядке
    sorted_regions_counter = dict(sorted(regions_counter.items(),
                                         key=lambda item: item[1],
                                         reverse=True
                                         )
                                  )

    # Выводим результаты
    # print("Unique regions count:")
    # for region, count in sorted_regions_counter.items():
    #     print(f"{region}: {count}")
    return sorted_regions_counter


def write_regions_to_file(regions_counter, output_file_path):
    with open(output_file_path, 'w', newline='') as output_file:
        csv_writer = csv.writer(output_file)
        for region, count in regions_counter.items():
            # Записываем одну строку для каждого региона
            row = list(region) + [count]
            csv_writer.writerow(row)


def read_txt_file(file_path):
    data_dict = {}
    with open(file_path, 'r') as file:
        # Пропускаем заголовок
        next(file)
        for line in file:
            # Разбиваем строку по пробелам
            parts = line.strip().split()
            # Получаем ключ (первая колонка) и значение (пятая колонка)
            key = parts[0]
            value = parts[3]
            # Добавляем в словарь
            data_dict[key] = value
        data_dict['SINGLETONE'] = 'singletone'
    return data_dict


def process_file(input_file_path, output_file_path, data_dict):
    with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
        for line in input_file:
            keys = line.strip().split(',')
            # Игнорируем последнее число в строке
            count = keys[-1]
            keys = keys[:-1]
            result_line = []
            for key in keys:
                if key in data_dict:
                    result_line.append(data_dict[key])
            # Записываем строку в выходной файл с добавлением последнего числа из словаря
            output_file.write(','.join(result_line) + ',' + count + '\n')


input_proteins_table = sys.argv[1]
small_cluster_size = sys.argv[2]
output_result_file_name = sys.argv[3]

update_second_column(input_proteins_table, 'output_step_1.txt')

main('output_step_1.txt',
     'output_step_2.csv',
     small_cluster_threshold=int(small_cluster_size))

with open('output_step_2.csv', 'r') as f_in, open('output_step_3.csv', 'w', newline='') as f_out:
    reader = csv.reader(f_in)
    writer = csv.writer(f_out)
    for line in reader:
        reversed_line = reverse_proteins(','.join(line))
        writer.writerow(reversed_line.split(','))

filter_csv('output_step_3.csv', 'output_step_4.csv')

# Использование функции для обработки CSV файла
process_csv('output_step_4.csv', 'output_step_5.csv')

regions_counter = count_unique_regions('output_step_5.csv')

write_regions_to_file(regions_counter, 'output_step_6.csv')

data_dict = read_txt_file(f'output_step_1.txt')

process_file('output_step_6.csv', output_result_file_name, data_dict)
