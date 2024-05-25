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
    # Find the second occurrence of the "_" character and extract the number after it
    second_underscore_index = protein_string.find('_', protein_string.find('_') + 1)
    number_string = protein_string[second_underscore_index + 1:]
    # Convert the string to a number, taking into account leading zeros
    return int(number_string.lstrip('0'))


def write_to_output_file(clust_region_list, output_file):
    with open(output_file, 'a') as f:
        f.write(','.join(clust_region_list) + '\n')


def main(input_file, output_file, small_cluster_threshold=0):
    # Sort the file by the second column
    with open(input_file, 'r') as f:
        lines = f.readlines()
        header = lines[0]
        data = sorted(lines[1:], key=lambda x: x.split(' ')[2])

    clust_region = []
    prev_protein_number = None

    for line in data:
        columns = line.strip().split(' ')
        # Check if there are four columns
        if len(columns) < 3:
            print("Error: Invalid number of columns per row:", line)
            continue
        protein = columns[2]
        size_clust = int(columns[1])
        clust_region_value = columns[0]
        if size_clust <= small_cluster_threshold:
            clust_region_value = 'x' + clust_region_value

        # Get the number from the second column (protein)
        try:
            protein_number = get_number_from_protein(protein)
        except ValueError:
            print("Error: Cannot extract number from string:", protein)
            continue

        if prev_protein_number is None:
            prev_protein_number = protein_number

        # Check the difference between the current and previous numbers
        if abs(protein_number - prev_protein_number) <= 1:
            # Add a value to the clust_region sheet
            # Take the value from the last column
            clust_region.append(clust_region_value + '_' + columns[-2])
        else:
            # The difference is greater than 1, write the current list to a new file
            write_to_output_file(clust_region, output_file)
            clust_region = []
            # Add the current value to a new list
            # Take the value from the last column
            clust_region.append(clust_region_value + '_' + columns[-2])

        prev_protein_number = protein_number
        
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
    # Remove the "_inner" and "_brex" annotations and sort the protein names
    processed_names = [name.replace("_inner", "").replace("_brex", "") for name in names]
    return tuple(processed_names)


def count_unique_regions(file_path):
    regions_counter = Counter()

    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            # Convert protein names and count unique regions
            region = process_protein_names(row)
            regions_counter[region] += 1

    # Sort the dictionary by values (number of occurrences of the region) in descending order
    sorted_regions_counter = dict(sorted(regions_counter.items(),
                                         key=lambda item: item[1],
                                         reverse=True
                                         )
                                  ))
    return sorted_regions_counter


def write_regions_to_file(regions_counter, output_file_path):
    with open(output_file_path, 'w', newline='') as output_file:
        csv_writer = csv.writer(output_file)
        for region, count in regions_counter.items():
            # Write one line for each region
            row = list(region) + [count]
            csv_writer.writerow(row)


def read_txt_file(file_path):
    data_dict = {}
    with open(file_path, 'r') as file:
        next(file)
        for line in file:
            parts = line.strip().split()
           # Get the key (first column) and value (fifth column)
            key = parts[0]
            value = parts[3]
            data_dict[key] = value
        data_dict['SINGLETONE'] = 'singletone'
    return data_dict


def process_file(input_file_path, output_file_path, data_dict):
    with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
        for line in input_file:
            keys = line.strip().split(',')
            # Ignore the last number in the line
            count = keys[-1]
            keys = keys[:-1]
            result_line = []
            for key in keys:
                if key in data_dict:
                    result_line.append(data_dict[key])
            # Write a string to the output file adding the last number from the dictionary
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

# Using a function to process a CSV file
process_csv('output_step_4.csv', 'output_step_5.csv')

regions_counter = count_unique_regions('output_step_5.csv')

write_regions_to_file(regions_counter, 'output_step_6.csv')

data_dict = read_txt_file(f'output_step_1.txt')

process_file('output_step_6.csv', output_result_file_name, data_dict)
