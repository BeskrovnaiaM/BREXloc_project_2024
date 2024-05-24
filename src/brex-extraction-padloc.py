import os
import pandas as pd
import sys

from Bio import SeqIO


def count_brex(file_path):
    """
    This function reads a CSV file at the given file path, extracts unique system numbers and system types
    containing 'brex', and returns a list of system numbers and their corresponding system types.

    Parameters:
    file_path (str): The file path to the CSV file.

    Returns:
    tuple: Two lists - system_num_list containing unique system numbers, and system_type_list containing
           corresponding system types with 'brex'.
    """
    df = pd.read_csv(file_path, header=None, skiprows=1, sep=',')

    unique_system_num_values = df[0].unique()
    unique_system_type_values = df[df[2].str.contains('brex')][2].unique()

    system_num_list = []
    system_type_list = []

    for val in unique_system_num_values:
        corresponding_val = df[(df[0] == val) & (df[2].str.contains('brex'))][2].values
        if len(corresponding_val) > 0:
            system_num_list.append(val)
            system_type_list.append(corresponding_val[0])

    return system_num_list, system_type_list


def process_csv_file(file_path, system_num, system_type):
    """
    This function processes a CSV file at the given file path based on the system number and system type provided.
    It extracts specific information such as the first and last numbers, corresponding names, and coordinates.

    Parameters:
    file_path (str): The file path to the CSV file.
    system_num (int): The system number to filter the data.
    system_type (str): The system type to match in the data.

    Returns:
    tuple: Contains first_number, last_number, name_contig_first, name_contig_last, start_coords, and last_coords.
    """
    regex = r"\b" + system_type + r"\b"
    system_num_str = str(system_num)
    df = pd.read_csv(file_path, header=None, sep=',')
    first_brex_index = df[(df[0] == system_num_str) & (df[2].str.match(regex))].index[0]
    last_brex_index = df[(df[0] == system_num_str) & (df[2].str.match(regex))].index[-1]
    first_number = df.iloc[first_brex_index, 11]
    last_number = df.iloc[last_brex_index, 12]
    name_contig_first = df.iloc[first_brex_index, 1]
    name_contig_last = df.iloc[last_brex_index, 1]
    matching_rows = df[(df[0] == system_num_str) & (df[2].str.match(regex))]
    start_coords = matching_rows[11].tolist()
    last_coords = matching_rows[12].tolist()
    return first_number, last_number, name_contig_first, name_contig_last, start_coords, last_coords


def find_closest_rows(file_path,
                      first_number,
                      last_number,
                      name_contig_first,
                      name_contig_last
                      ):
    """
    This function finds the closest rows in a file based on specified criteria.

    Parameters:
    file_path (str): The file path to the file.
    first_number (int): The first number for comparison.
    last_number (int): The last number for comparison.
    name_contig_first (str): The name for the first comparison.
    name_contig_last (str): The name for the last comparison.

    Returns:
    tuple: Contains closest_first, closest_last, coord_first, and coord_last.
    """
    first_number = int(first_number)
    last_number = int(last_number)
    first_number -= 10000
    last_number += 10000
    closest_first = None
    closest_last = None
    min_diff_first = float('inf')
    min_diff_last = float('inf')
    if first_number - 10000 < 0:
        first_number = 0

    with open(file_path, 'r') as file:
        line_number = 0
        for line in file:
            cols = line.strip().split('\t')
            try:
                col4 = int(cols[3])
                col5 = int(cols[4])
                col1 = cols[0]
                diff_first = abs(first_number - col4)
                diff_last = abs(last_number - col5)

                if col4 >= first_number and diff_first < min_diff_first and col1 == name_contig_first:
                    closest_first = line_number
                    min_diff_first = diff_first
                    coord_first = col4

                if col5 <= last_number and diff_last < min_diff_last and col1 == name_contig_last:
                    closest_last = line_number
                    min_diff_last = diff_last
                    coord_last = col5
            except ValueError:
                pass

            line_number += 1

    return closest_first, closest_last, coord_first, coord_last


def write_flank_gff(file_path_to_gff,
                    name_assemb,
                    name_contig_first,
                    num,
                    ident_type,
                    closest_first,
                    closest_last,
                    output_dir
                    ):
    """
    This function writes specific lines from a GFF file to a new GFF file based on the provided parameters.

    Parameters:
    file_path_to_gff (str): The file path to the GFF file.
    name_assemb (str): The name of the assembly.
    name_contig_first (str): The name of the first contig.
    num (int): A number for identification.
    ident_type (str): The type for identification.
    closest_first (int): The index of the closest first line.
    closest_last (int): The index of the closest last line.

    Returns:
    str: The file path to the newly created GFF file containing the selected lines.
    """

    name_output_gff = (name_assemb + '_' +
                       ident_type + '_' +
                       str(num) + '_' +
                       name_contig_first +
                       '_reg.gff'
                       )

    output_gff = os.path.join(output_dir, name_output_gff)

    closest_last += 1
    lines_to_read = range(closest_first, closest_last)
    result = [line for i, line in enumerate(open(file_path_to_gff)) if i in lines_to_read]

    with open(output_gff, 'a') as file:
        for string in result:
            file.write(string)
    return output_gff


def check_genes_positions(input_file_path,
                          first_gene,
                          last_gene,
                          start_coords,
                          last_coords
                          ):
    """
    This function checks gene positions in a GFF file based on specified criteria and updates a new column accordingly.

    Parameters:
    input_file_path (str): The file path to the input GFF file.
    first_gene (int): The position of the first gene.
    last_gene (int): The position of the last gene.
    start_coords (list): List of start coordinates.
    last_coords (list): List of last coordinates.

    Returns:
    None
    """
    df = pd.read_csv(input_file_path, header=None, sep='\t')
    df['new_col'] = 'inner'
    start_coords = [int(coord) for coord in start_coords]
    last_coords = [int(coord) for coord in last_coords]
    strand_value = ''
    for index, row in df.iterrows():
        if row[3] in start_coords and row[4] in last_coords:
            df.at[index, 'new_col'] = 'brex'
            strand_value = row[6]
    for index, row in df.iterrows():
        if row[3] >= int(first_gene) and row[4] <= int(last_gene) and not row['new_col'] == 'brex':
            df.at[index, 'new_col'] = 'inner'
        elif row[3] < int(first_gene) and row[4] < int(first_gene):
            if strand_value == '+':
                df.at[index, 'new_col'] = 'upstream'
            else:
                df.at[index, 'new_col'] = 'downstream'
        elif row[3] > int(last_gene) and row[4] > int(last_gene):
            if strand_value == '+':
                df.at[index, 'new_col'] = 'downstream'
            else:
                df.at[index, 'new_col'] = 'upstream'

    output_file_path = input_file_path.replace('.gff', '_positions.gff')
    df.to_csv(output_file_path, sep='\t', header=False, index=False)


def extract_sequences(input_file,
                      output_file,
                      start_line,
                      end_line,
                      name_assemb,
                      name_contig_first,
                      num,
                      ident_type,
                      output_dir
                      ):
    """
    Extracts sequences from a FASTA file based on the specified start and end lines.

    Parameters:
    input_file (str): The input FASTA file path.
    output_file (str): The output file path to write the extracted sequences.
    start_line (int): The start line for extraction.
    end_line (int): The end line for extraction.
    name_assemb (str): The name of the assembly.
    name_contig_first (str): The name of the first contig.
    num (int): A number for identification.
    ident_type (str): The type for identification.

    Returns:
    None
    """
    sequences = []

    name_output_faa = (name_assemb + '_' +
                       ident_type + '_' +
                       str(num) + '_' +
                       name_contig_first +
                       '_reg.faa'
                       )

    output_faa = os.path.join(output_dir, name_output_faa)

    with open(input_file, 'r') as file:
        for i, record in enumerate(SeqIO.parse(file, 'fasta')):
            if start_line <= i <= end_line:
                sequences.append(record)

    with open(output_file, 'a') as out_file:
        SeqIO.write(sequences, out_file, 'fasta')

    with open(output_faa, 'a') as file:
        SeqIO.write(sequences, file, 'fasta')


def extract_coord_loc(name_assemb,
                      num,
                      ident_type,
                      name_contig_first,
                      name_contig_last,
                      first_gene,
                      last_gene,
                      closest_first,
                      closest_last,
                      coord_first,
                      coord_last,
                      coor_table
                      ):
    """
    Writes coordinate location information to a coordinate table.

    Parameters:
    name_assemb (str): The name of the assembly.
    num (int): A number for identification.
    ident_type (str): The type for identification.
    name_contig_first (str): The name of the first contig.
    name_contig_last (str): The name of the last contig.
    first_gene (int): The position of the first gene.
    last_gene (int): The position of the last gene.
    closest_first (int): The index of the closest first line.
    closest_last (int): The index of the closest last line.
    coord_first (int): The coordinate of the first location.
    coord_last (int): The coordinate of the last location.
    coor_table (str): The file path to the coordinate table.

    Returns:
    None
    """
    with open(coor_table, 'a') as file:
        file.write(f'{name_assemb}\t{ident_type}\t{num}\t'
                   f'{name_contig_first}\t{name_contig_last}\t'
                   f'{first_gene}\t{last_gene}\t{closest_first}\t'
                   f'{closest_last}\t{coord_first}\t{coord_last}\n'
                   )


def extract_faa(input_dir_path,
                sample_names,
                output_file_faa,
                coor_table,
                output_dir
                ):
    """
    Extracts information and sequences from multiple files in a directory and writes them to output files.

    Parameters:
    input_dir_path (str): The path to the input directory.
    sample_names (list): List of sample names.
    output_file_faa (str): The output file path for writing the extracted sequences.
    coor_table (str): The output file path to the coordinate table.

    Returns:
    list: List of files with errors encountered during processing.
    """
    file_num = 1
    errors_csv = []
    for file in sample_names:
        file_path_to_csv = os.path.join(input_dir_path, file + '_padloc.csv')
        name_assemb = file[0:15]
        print(file_path_to_csv)
        system_num_list, system_type_list = count_brex(file_path_to_csv)
        if len(system_type_list) == 0:
            errors_csv.append(file)
            print(f'Oops, there is no brex in file {file_num}')
        else:
            for num, sys_type in zip(system_num_list, system_type_list):
                (first_gene, last_gene,
                 name_contig_first, name_contig_last,
                 start_coords, last_coords) = process_csv_file(file_path_to_csv, num, sys_type)
                if name_contig_first != name_contig_last:
                    errors_csv.append(file)
                    print(f'Brex system {file_num} is broken. Two contigs')
                    break
                file_path_to_gff = os.path.join(input_dir_path, file + '_prodigal.gff')
                (closest_first, closest_last,
                 coord_first, coord_last) = find_closest_rows(file_path_to_gff,
                                                              first_gene,
                                                              last_gene,
                                                              name_contig_first,
                                                              name_contig_last
                                                              )
                file_path_to_faa = os.path.join(input_dir_path, file + '_prodigal.faa')
                input_file_path_to_flank_gff = write_flank_gff(file_path_to_gff,
                                                               name_assemb,
                                                               name_contig_first,
                                                               num,
                                                               sys_type,
                                                               closest_first,
                                                               closest_last,
                                                               output_dir
                                                               )
                check_genes_positions(input_file_path_to_flank_gff,
                                      first_gene,
                                      last_gene,
                                      start_coords,
                                      last_coords
                                      )
                extract_sequences(file_path_to_faa,
                                  output_file_faa,
                                  closest_first,
                                  closest_last,
                                  name_assemb,
                                  name_contig_first,
                                  num,
                                  sys_type,
                                  output_dir
                                  )
                extract_coord_loc(name_assemb,
                                  num,
                                  sys_type,
                                  name_contig_first,
                                  name_contig_last,
                                  first_gene,
                                  last_gene,
                                  closest_first,
                                  closest_last,
                                  coord_first,
                                  coord_last,
                                  coor_table
                                  )
                print(f'I am done with file number {file_num}')
        file_num += 1
    return errors_csv


input_padloc_dir = sys.argv[1]
output_results_dir_path = sys.argv[2]
output_fasta = sys.argv[3]
output_coor_table = sys.argv[4]

if not os.path.exists(output_results_dir_path):
    os.mkdir(output_results_dir_path)

pattern = 'prodigal.gff'

samples = set()

for paths, dirs, files in os.walk(input_padloc_dir):
    for file in files:
        if pattern in file:
            sample = file.split('_prodigal.gff')[0]
            samples.add(sample)

genomes_wo_brex = extract_faa(input_padloc_dir,
                              samples,
                              output_fasta,
                              output_coor_table,
                              output_results_dir_path)

for i in genomes_wo_brex:
    with open("genomes_wo_brex.txt", 'a') as file:
        file.write(f'{i}\n')

