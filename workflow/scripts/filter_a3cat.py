import pandas as pd
import numpy as np


# Read the input files
big_table_df = pd.read_csv(snakemake.input['a3cat_table'], sep='\t', low_memory=False)

with open(snakemake.input['proteome_list'], 'r') as file:
    # Extract genome ID pattern
    genome_ids = ['_'.join(line.strip().split('_')[:2]) for line in file]

# Filter rows where Genbank Accession matches genome_id in the proteome_list file
filtered_df = big_table_df[big_table_df['Genbank Accession'].isin(genome_ids) | big_table_df['Refseq Accession'].isin(genome_ids)]


# Create "RefSeq/GenBank id" column which take RefSeq Id if available, and Genbank otherwise
filtered_df = filtered_df.copy()
filtered_df['RefSeq/GenBank Accession'] = np.where((filtered_df['Refseq Accession'].notnull()) & (filtered_df['Refseq Accession'] != ''), filtered_df['Refseq Accession'], filtered_df['Genbank Accession'])


# Select the specified columns
columns_to_select = [
    'RefSeq/GenBank Accession', 'Submission Date', 'TaxId', 'SubPhylum', 'Order', "Genus", "Species",
    'Total length', 'Number of scaffolds', 'Scaffold L50', 'Scaffold N50', "Annotation Date",
    'Total genes', 'Arthropoda Complete BUSCO', 'Arthropoda Single BUSCO', 'Arthropoda Duplicated BUSCO',
    'Arthropoda Fragmented BUSCO', 'Arthropoda Missing BUSCO'
]


#columns_to_select = [
#    'Genbank Accession', 'Refseq Accession', 'Submission Date', 'Bioproject', 'Organism Name', 'Strain', 'TaxId', 'SubPhylum', 'Order', 'Assembly Level',
#    'Total length', 'Number of contigs', 'Contig N50', 'Contig L50', 'Number of scaffolds',
#    'Scaffold N50', 'Scaffold L50', 'Total genes', 'Protein-coding genes', 'Pseudogenes', 'Non-coding genes', 'Arthropoda Complete BUSCO', 
#    'Arthropoda Single BUSCO', 'Arthropoda Duplicated BUSCO', 'Arthropoda Fragmented BUSCO',
#    'Arthropoda Missing BUSCO'
#]
selected_data = filtered_df[columns_to_select]

# Write the output file
selected_data.to_csv(snakemake.output['selected_data'], sep='\t', index=False)

