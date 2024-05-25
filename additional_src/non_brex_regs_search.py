"""
    Counts the total number of non-BREX systems
    taking into account them localization

    Input positional args:
    path to directory with PADLOC csv-files
    path to directory with gff-files with localization marks
    
    Output files:
    'Total_non_brex_systems.tsv'
    'Total_inner_systems.tsv'
    'Total_other_systems.tsv'
    'Top_10_inner_systems.png'
"""

import os
import sys
import pandas as pd


def write_results(output_file_name: str, resuls_to_write: dict) -> None:
    """
    Takes as input a dictionary with counted ES and writes it to a file
    """
    with open(output_file_name, mode='w') as file:
        file.write('System_type' + '\t' + 'Number' + '\n')
        for key, value in resuls_to_write.items():
            new_line = key + '\t' + str(value)  +'\n'
            file.write(new_line)

def plot_sys(syst_dict: dict, n_to_choose: int = 10):
    """Selects the specified top number of defense systems 
    from the dictionary in the form of a dataframe"""
    
    df_plot = pd.DataFrame.from_dict(syst_dict, orient='index')
    df_plot.rename(columns={0: 'Number'}, inplace=True)
    top_n = df_plot.sort_values(by=['Number'], ascending=False).iloc[:n_to_choose]
    return top_n   

# input_padloc_dir = sys.argv[1]
# input_pos_dir = sys.argv[2]

input_padloc_dir = '/home/holydiver/Main/IB/Prrr/Data/03_Padloc_reannotation/Reann'
input_pos_dir = '/home/holydiver/Main/IB/Prrr/Data/02_BREX_regs/3_BREX_regs'

total_non_brex_def_systems = {}
total_inner_def_systems = {}
total_other_def_systems = {}

loci_with_inner_sys = 0

gff_header = ['seqid', '_2', '_3', '_4', '_5', '_6', '_7', '_8', 'Comment', 'Loc']

for file_name in os.listdir(input_padloc_dir):

    curr_non_brex_def_systems = {}
    curr_inner_def_systems = {}
    curr_other_def_systems = {}

    if file_name.endswith('csv'):
        df = pd.read_csv(os.path.join(input_padloc_dir, file_name)).iloc[:, :4]

        file_name_pos = file_name.split('padloc')[0] + 'positions.gff'
        df_pos = pd.read_csv(os.path.join(input_pos_dir, file_name_pos),
                             sep='\t', names=gff_header).iloc[:, [0, 8, 9]]

        # creating a unique ID allows you to take into account matching names of system types
        df['Name'] = df['system'] + '%' + df['system.number'].astype('str')

        df_pos['target.name'] = (df_pos['seqid'] +
                                 '_' +
                                 [i for i in map(lambda x: x.split(';')[0].split('_')[1], df_pos['Comment'])])
        
        df_non_brex = df[df['system'].str[:4] != 'brex'].merge(df_pos.iloc[:, [2, 3]], how='left')
        uniq_sys = df_non_brex.loc[:, ('Name', 'Loc')].groupby('Name').agg(lambda x: set(x))

        # Results lists
        curr_non_brex_def_systems = [i for i in map(lambda x: x.split('%')[0], uniq_sys.index.unique())]
        curr_inner_def_systems = [i for i in
                                  map(lambda x: x.split('%')[0], uniq_sys[uniq_sys['Loc'] == {'inner'}].index.unique())]
        curr_other_def_systems = [i for i in
                                  map(lambda x: x.split('%')[0], uniq_sys[uniq_sys['Loc'] != {'inner'}].index.unique())]
                
        if len(curr_inner_def_systems) > 0:
            loci_with_inner_sys += 1
                    
        # Counting
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

#  Outputs
output_files_names = ('Total_non_brex_systems.tsv', 'Total_inner_systems.tsv', 'Total_other_systems.tsv')
def_systems = (total_non_brex_def_systems, total_inner_def_systems, total_other_def_systems)

for file_name, sys_dict in zip(output_files_names, def_systems):
    write_results(file_name, sys_dict)

