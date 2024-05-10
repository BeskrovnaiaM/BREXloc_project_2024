import os
import re
import sys

from Bio import SeqIO


def find_cds_in_range(gbk_file, protocore_list):
    cds_list = []
    brex_names = ['BrxA', 'BrxB', 'BrxC', 'BrxD', 'BrxE', 'BrxF', 'BrxHI',
                  'BrxHII', 'BrxL', 'BrxP', 'PglW', 'PglX', 'PglXI', 'PglZ']
    pattern = re.compile(r'\[(\d+):(\d+)\]\(([\+\-])\)')

    for core in protocore_list:
        contig = core['Contig']
        brex_type = core['Product']
        start = core['Start']
        end = core['End']
        for record in SeqIO.parse(gbk_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    if "locus_tag" in feature.qualifiers:
                        locus_tag = feature.qualifiers["locus_tag"][0]
                    else:
                        locus_tag = "Unknown"
                    if "gene_functions" in feature.qualifiers:
                        gene_name = feature.qualifiers["gene_functions"][0]
                    else:
                        gene_name = "Unknown"
                    match = re.match(pattern, str(feature.location))
                    if match:
                        start_cds = int(match.group(1))
                        end_cds = int(match.group(2))
                        strand = match.group(3)
                    if start_cds >= start and end_cds <= end:
                        if any(name in gene_name for name in brex_names):
                            protein_type = 'brex'
                        else:
                            protein_type = 'inner'
                        cds_list.append({"Contig": contig,
                                         "Brex_type": brex_type,
                                         "Start_reg": start,
                                         "End_reg": end,
                                         "Locus_tag": locus_tag,
                                         "Product": gene_name,
                                         "Protein": protein_type,
                                         "Start": start_cds,
                                         "End": end_cds,
                                         "Strand": strand,
                                         "Translation": feature.qualifiers['translation']
                                         })
    return cds_list


def find_brex_protocore(gbk_file):
    protocore_list = []
    pattern = re.compile(r'\[(\d+):(\d+)\]\(([\+\-])\)')
    contig_region_name = os.path.basename(gbk_file)[:-4]
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type == "proto_core":
                product = feature.qualifiers['product'][0]
                if 'brex' in product:
                    match = re.match(pattern, str(feature.location))
                    if match:
                        start_core = int(match.group(1))
                        end_core = int(match.group(2))
                    protocore_list.append({"Contig": contig_region_name,
                                           "Product": product,
                                           "Start": start_core,
                                           "End": end_core
                                           })
    return protocore_list


def write_gff_and_fasta(cds_list, folder_path, fasta_filename="output.fasta"):
    contig = cds_list[0]['Contig']

    if not os.path.exists(os.path.join(folder_path, 'regions_gff')):
        os.mkdir(os.path.join(folder_path, 'regions_gff'))
    gff_filename = contig + ".gff"
    gff_path = os.path.join(folder_path, 'regions_gff', gff_filename)

    # Write GFF file
    with open(gff_path, 'w') as gff_file:
        header = "\t".join(['Contig',
                            'Brex_type',
                            'Start_reg',
                            'End_reg',
                            'Locus_tag',
                            'Product',
                            'Start',
                            'End',
                            'Strand',
                            'Protein'
                            ])
    # Запись заголовка в файл
        gff_file.write(header + '\n')
        for entry in cds_list:
            gff_line = "\t".join([entry['Contig'],
                                  entry['Brex_type'],
                                  str(entry['Start_reg']),
                                  str(entry['End_reg']),
                                  entry['Locus_tag'],
                                  entry['Product'],
                                  str(entry['Start']),
                                  str(entry['End']),
                                  entry['Strand'],
                                  entry['Protein']
                                  ])
            gff_file.write(gff_line + '\n')

    # Write or append to FASTA file
    with open(os.path.join(folder_path, fasta_filename), 'a') as fasta_file:
        for entry in cds_list:
            protein_name = (entry['Contig'] + '_' +
                            entry['Locus_tag'] + '_' +
                            entry['Brex_type'] + '_' +
                            entry['Protein'])
            fasta_seq = entry['Translation'][0]
            fasta_entry = SeqIO.SeqRecord(fasta_seq, id=protein_name, description='')
            SeqIO.write(fasta_entry, fasta_file, 'fasta')


def get_file_names(directory):
    file_names = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_names.append(os.path.join(root, file))
    return file_names


def gbk_parser(input_files, folder_path, output_fasta='output.fasta'):
    file_names = get_file_names(input_files)
    
    # Удаление файлов, не оканчивающихся на .gbk
    file_names = [file for file in file_names if file.endswith('.gbk')]
    
    # Получение полного пути к файлу, который нужно исключить
    input_folder_name = os.path.basename(os.path.normpath(input_files))
    forbidden_file_path = os.path.join(input_files, input_folder_name + '.gbk')
    
    # Удаление файла, полный путь которого совпадает с forbidden_file_path
    file_names = [file for file in file_names if file != forbidden_file_path]

    for file in file_names:
        protocore_list = find_brex_protocore(file)
        if protocore_list:
            cds_list = find_cds_in_range(file, protocore_list)
            write_gff_and_fasta(cds_list, folder_path, output_fasta)


input_gbk_dir_path = sys.argv[1]
output_results_dir_path = sys.argv[2]

if __name__ == "__main__":
    gbk_parser(input_gbk_dir_path, output_results_dir_path)
