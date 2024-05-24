import os
import pandas as pd
import re
import sys


def select_gene_name(comment_line: str, pattern: str) -> str:
    """Search gene name in gff-file and return it"""

    gene_name = re.search(pattern, comment_line).group(2)
    return gene_name


def choose_annotation(df: pd.DataFrame,
                      annot_heterog_fraction: float = 0.3,
                      no_annot_fractions: float = 0.5) -> str:
    """
    Choose which annotation (Padloc or HH_Suite) be set as final.
    Padloc is preferential as more specific to defense systems

    Choose HH-Suite when:
    - Any Padloc annotation assign to no more than 'annot_heterog_fraction'
      (high heterogeneity of annotation)
    - More than 'no_annot_fractions' of the cluster have missing Padloc annotation
    """

    check_if_no_annotation = df.Padloc_ann.value_counts().index[0] == 'no_ann'
    max_fraction = (df.Padloc_ann.value_counts() / df.shape[0]).max()

    if max_fraction <= annot_heterog_fraction:
        annot = df.HHSuite_ann.values[0]
    elif check_if_no_annotation:
        if max_fraction > no_annot_fractions:
            annot = df.HHSuite_ann.values[0]
        else:
            second_fraction = (df.Padloc_ann.value_counts() / df.shape[0]).iloc[1]
            if second_fraction <= annot_heterog_fraction:
                annot = df.HHSuite_ann.values[0]
            else:
                annot = df.Padloc_ann.value_counts().index[1]
    else:
        annot = df.Padloc_ann.value_counts().index[0]

    return annot


# Args
input_clusters_list = sys.argv[1]
input_padloc_dir = sys.argv[2]
input_gff_dir = sys.argv[3]
input_clusters_size = sys.argv[4]
input_hh_suite_annotations = sys.argv[5]
output_tables_prefix = sys.argv[6]

# Main table
df_clusters = pd.read_csv(input_clusters_list, 
                          sep='\t', 
                          names=['Cluster', 'Protein']
                          )

# Add Padloc annotations
df_clusters.set_index('Protein', inplace=True)
df_clusters['Padloc_ann'] = 'no_ann'

for file_name in os.listdir(input_padloc_dir):
    if file_name.endswith('csv'):
        df_annot = pd.read_csv(os.path.join(input_padloc_dir, 
                                            file_name
                                            )
                               ).iloc[:, :12]
        
        df_annot_sel = (df_annot.groupby('target.name')
                                .agg('min')
                                .loc[:, 'protein.name'])

        idxs_intersection = list(set(df_annot_sel.index) & set(df_clusters.index))

        df_clusters.loc[idxs_intersection, 'Padloc_ann'] = df_annot_sel.loc[idxs_intersection].values

# Add localisation marks and brex-system types
gff_cols_names = ('Chrom', '_2', '_3', '_4', '_5', '_6', '_7', '_8', 'Comment', 'Loc')
pattern_for_gene = re.compile(r'(ID=\d+_(\d+))')
df_clusters['Localisation'] = 'no_loc'
df_clusters['Brex_type'] = 'no_type'

for file_name in os.listdir(input_gff_dir):
    if file_name.endswith('positions.gff'):
        df_loc = pd.read_csv(os.path.join(input_gff_dir, file_name),
                             sep='\t',
                             names=gff_cols_names)

        brex_locus_type = '_'.join(file_name.split('_')[2:5])

        df_loc['Protein'] = (df_loc['Chrom'] +
                             '_' +
                             df_loc['Comment'].apply(select_gene_name,
                                                     args=(pattern_for_gene,)
                                                     )
                             )

        df_loc.set_index('Protein', inplace=True)

        idxs_intersection = list(set(df_loc.index) & set(df_clusters.index))

        df_clusters.loc[idxs_intersection, 'Localisation'] = df_loc['Loc'].loc[idxs_intersection].values
        df_clusters.loc[idxs_intersection, 'Brex_type'] = brex_locus_type

# Add cluster sizes
df_clusters_size = pd.read_csv(input_clusters_size, sep='\t', names=['Cluster', 'Cluster_size'])

df_clusters['Protein'] = df_clusters.index

df_clusters = pd.merge(df_clusters, df_clusters_size, how='left', on='Cluster')

# Add HH-Suite annotations
df_clusters['HHSuite_ann'] = 'no_ann'

pattern_ann = re.compile(r'(\S+) ; (\S+) ;')

only_first_annotations = {}  # Top-1 annotations for main table

for file_name in os.listdir(input_hh_suite_annotations):
    if file_name.endswith('hhr'):
        with open(os.path.join(input_hh_suite_annotations, file_name)) as file:
            cluster_id = file_name[:-4]

            count = 0
            while count < 9:
                file.readline()
                count += 1

            curr_line = file.readline()
            curr_first_annotation = re.search(pattern_ann, curr_line).group(2)
            only_first_annotations[cluster_id] = curr_first_annotation

only_first_annotations = pd.Series(only_first_annotations)
df_clusters.set_index('Cluster', inplace=True)
df_clusters.loc[only_first_annotations.index, 'HHSuite_ann'] = only_first_annotations

# Choose one annotation per cluster
df_annotation = pd.DataFrame(df_clusters.groupby('Cluster')[['Padloc_ann', 'HHSuite_ann']]
                                        .apply(choose_annotation),
                             columns=['Annotation']
                             ).reset_index()

df_clusters = pd.merge(df_clusters,
                       df_annotation,
                       how='left',
                       on='Cluster'
                       )

# Writing final tables
output_all_annot_name = output_tables_prefix + '_all_ann_clusters_table.tsv'
output_one_annot_name = output_tables_prefix + '_final_clusters_table.tsv'

all_annotations_colnames = ['Cluster', 'Cluster_size', 'Protein', 'Padloc_ann',
                            'HHSuite_ann', 'Annotation', 'Localisation', 'Brex_type'
                            ]

one_annotation_colnames = ['Cluster', 'Cluster_size', 'Protein',
                           'Annotation', 'Localisation', 'Brex_type'
                           ]

df_clusters.to_csv(output_all_annot_name,
                   sep='\t',
                   header=True,
                   index=False,
                   columns=all_annotations_colnames
                   )

df_clusters.to_csv(output_one_annot_name,
                   sep='\t',
                   header=True,
                   index=False,
                   columns=one_annotation_colnames
                   )

