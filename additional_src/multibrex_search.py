"""
Checks the uniqueness of brex systems using PADLOC csv-files.
Saves a list of all multi-brex files to the table 'multi_brex_files.tsv'

Can create a copy/overwrite PADLOC csv-files,
adding a column with a unique ID: system_system.number_seqid  
  
  
Args:
1 - list of names of PADLOC csv-files
It is assumed that all files are in the same directory
File names are written to a summary table with a list of multi-brex files
2 - path to the directory with these files
3 optional - specify the folder name if needed to save updated files with an augmented column
If  specify REWRITE, it will add a column, overwriting the original csv-files
"""

import os
import sys
import pandas as pd


def myf(file_path, file_name, write_csv=False):
    """
    From the 'system' column selects only 'brex'.
    Counts how many unique 'system.number' there are.
    If greater than 1, returns 
    the file name, 
    the number of brex systems,
    number of contigs for brex, 
    and a list of brex systems with names like: system_system.number_seqid
    
    Saves a csv file with a 'Name' column at the end with the names,
    if the 'write_csv' specified
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


input_files_list = sys.argv[1]
input_padloc_dir = sys.argv[2]

if len(sys.argv) == 4:
    if sys.argv[3] == 'REWRITE':
        output_dir = input_padloc_dir
    else:
        output_dir = sys.argv[3]
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
else:
    output_dir = False

files_list = []

with open(input_files_list) as file:
    for line in file:
        if not line.strip() == '':
            files_list.append(line.strip())

list_multi_brex_files = []

for name in files_list:
    res = myf(os.path.join(input_padloc_dir, name), name, output_dir)
    if res:
        list_multi_brex_files.append(res)

# Outputs
header = 'Accession' + '\t' + 'Total_brex' + '\t' + 'Contigs_brex' + '\t' + 'Brex_list' + '\n'
with open(os.path.join('multi_brex_files.tsv'), mode='w') as file:
    file.write(header)
    for record in list_multi_brex_files:
        f_name, b_num, c_num, b_names = '.'.join(record[0].split('.')[:-1]), record[1], record[2], ','.join(record[3])
        new_rec = f_name + '\t' + b_num + '\t' + c_num + '\t' + b_names + '\n'
        file.write(new_rec)

