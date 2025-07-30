# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


#conda activate env_SeqIO

from Bio import SeqIO
import re
import pandas as pd
from tqdm import tqdm
import time
from Bio import Entrez
import urllib.error

#fasta_file = "/Users/sagane/Desktop/test/GCF_021234035.1_SC_F0-13Bv2_filtered.faa"
#fasta_file = "/Users/sagane/Desktop/test/GCA_004104545.1_Arma_vul_BF2787_filtered.faa"
#output_file = "/Users/sagane/Desktop/test/test2.tsv"

def process_result(result):
    """
    Extracts relevant metadata from edirect result (see function bellow). + retrieve prot len using a different query
    """
    # Regular expressions to extract specific information from the result string
    base_transcript_id_match = re.search(r'>lcl\|(.*?)_cds_', result)
    locus_tag_search = re.search(r'\[locus_tag=(.*?)\]', result)
    gene_name_search = re.search(r'\[gene=(.*?)\]', result)
    gene_id_search = re.search(r'db_xref=.*?GeneID:(\d+)', result)
    location_search = re.search(r'\[location=(.*?)\]', result)
    description_search = re.search(r'\[protein=(.*?)\]', result)
    protein_id_search = re.search(r'\[protein_id=(.*?)\]', result)

    # Retrieve  metadata
    gene_id = gene_id_search.group(1) if gene_id_search and gene_id_search.group(1) else 'NA'
    gene_name = gene_name_search.group(1) if gene_name_search else 'NA'
    protein_id = protein_id_search.group(1) if protein_id_search else 'NA'
    protein_description = description_search.group(1) if description_search else 'NA'
    locus_tag = locus_tag_search.group(1) if locus_tag_search else 'NA'

    # Retrieve transcript ID: Use locus_id and position for Genbank (lacks standard transcript IDs); use standard ID for others..
    if locus_tag_search and location_search:
        location = (location_search.group(1).replace('join(', '').replace('complement(', '').replace(')', '').replace('..', ':'))
        transcript_id = 'NaN'
        transcript_url_suffix = f"{base_transcript_id_match.group(1) if base_transcript_id_match else 'NA'}?location={location}&report=gbwithparts"

    else:
        transcript_id = base_transcript_id_match.group(1) if base_transcript_id_match else 'NA'
        transcript_url_suffix = f"{base_transcript_id_match.group(1) if base_transcript_id_match else 'NA'}?report=gbwithparts"

    # Extract protein length from the output (len of nucleotides divided by 3 )
    nucleotide_seq_match = re.search(r'\n(.*?)(?:\n>|$)', result, re.DOTALL)
    nucleotide_sequence = nucleotide_seq_match.group(1).replace('\n', '') if nucleotide_seq_match else ''
    prot_len = (len(nucleotide_sequence) // 3) - 1 if nucleotide_sequence else 'NA'
    # return everything:
    return {
        'protein_id': protein_id,
        'prot_len': prot_len,
        'transcript_id': transcript_id,
        'transcript_url_suffix': transcript_url_suffix,
        'gene_id': gene_id,
        'gene_name': gene_name,
        'locus_tag': locus_tag,
        'protein_description': protein_description
    }


def determine_origin_and_genome_id(fasta_file):
    """
    Determines the origin and extracts the genome ID from the fasta file name.
    """
    # Extract genome ID
    genome_id_match = re.search(r'(GCF|GCA)_[0-9]+\.[0-9]+', fasta_file)
    genome_id = genome_id_match.group() if genome_id_match else "Unknown"

    # Determine origin based on the genome ID prefix
    if genome_id.startswith('GCF'):
        origin = "NCBI RefSeq assembly"
    elif genome_id.startswith('GCA'):
        origin = "Submitted GenBank assembly"
    else:
        origin = "Unknown"

    return origin, genome_id



def process_fasta(fasta_file, output_file):
    """
    Processes each record in a FASTA file, retrieves metadata using Entrez Direct, and compiles the results into a DataFrame.
    """
    # Determine the origin( Genbank or Refseq) and genome id
    origin, genome_id = determine_origin_and_genome_id(fasta_file)
    Entrez.email = "sagane.joye@unil.ch"
    Entrez.api_key = "3fb72e2b36e135fabd41ca74a88a3bdd3109"
    Entrez.max_tries = 5  # Maximum number of retries (default is 3)
    Entrez.sleep_between_tries = 3  # Seconds to sleep between retries (default is 15)
    df = pd.DataFrame(columns=['protein_id', 'prot_len', 'transcript_id', 'transcript_url_suffix', 'gene_id', 'gene_name', 'locus_tag', 'protein_description'])
    i = 0
    # Process each record in the FASTA file
    with open(fasta_file, 'r') as file:
        for record in tqdm(SeqIO.parse(file, 'fasta'), desc="Processing FASTA"):
            protein_id = record.id
            success = False
            attempts = 0
            while not success and attempts < 5:  # Retry up to 5 times
                try:
                    # Retrieve metadata using Entrez Direct
                    handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta_cds_na", retmode="text")
                    result = handle.read()
                    handle.close()
                    success = True  # If this line is reached, efetch was successful
                except urllib.error.HTTPError as e:
                    if e.code == 429:  # Too Many Requests
                        print("Received HTTP 429 Error, waiting longer...")
                        time.sleep(3)  # Wait longer if hit by rate limit
                    else:
                        raise e  # Re-raise other HTTP errors
                except Exception as e:
                    print(f"Error fetching data for {protein_id}: {e}")
                    break  # Break out of the loop on other errors
                finally:
                    attempts += 1

            if not success:
                continue  # Skip this record if all attempts fail

            # Process the result and append it to the DataFrame
            row_data = process_result(result)
            row_data['origin'] = origin
            row_data['genome_id'] = genome_id
            #print(row_data)
            df = pd.concat([df, pd.DataFrame([row_data])], ignore_index=True)
            #i += 1
            #if i > 50:
            #    break
            time.sleep(0.25)  # Wait for 0.5 seconds before processing the next record
    df.to_csv(output_file, sep="\t", index=False)

# usage:
process_fasta(snakemake.input['filtered'], snakemake.output['metadata_table'])


