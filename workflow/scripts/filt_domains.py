import pandas as pd

def filter_domains(input_files, output_file):
    # Define column names
    column_names = ['Prot_IDs', 'MD5 digest', 'seq_len', 'analysis', 'accession', 
                    'desc_1', 'Start', 'End', 'Score', 'Status', 'Date', 
                    'IPR_ID', 'description', 'GO', 'Pathways']

    df = pd.read_csv(input_files, sep='\t', header=None, names=column_names)
    filtered_df = df[df["Status"].eq('T') & df['IPR_ID'].str.contains(r'IPR\d+', na=False)]
    final_df = filtered_df[['Prot_IDs', 'IPR_ID', 'Start', 'End', 'description']]
    final_df.to_csv(output_file, sep='\t', index=False, header=True)


filter_domains(snakemake.input['domains'], snakemake.output['filtered_domains'])




#import pandas as pd


# read interproscan output 
#column_names = ['Prot_IDs', 'MD5 digest', 'seq_len', 'analysis', 'accession', 'desc_1', 'Start', 'End', 'Score', 'Status', 'Date', 'Prot_IDs', 'description', 'GO', 'Pathways']
#df = pd.read_csv(snakemake.output['domains'], sep='\t', header=None, names=column_names)

# Filter rows where column has 'T' (for True ) and contains 'IPR' identifiers
#filtered_df = df[df["Status"].eq('T') & df['Prot_IDs'].str.contains(r'IPR\d+')]
#final_df = filtered_df[['Prot_IDs', 'IPR_ID', 'Start', 'End', 'description']]
#final_df.to_csv(, sep='\t', index=False, header=True)





