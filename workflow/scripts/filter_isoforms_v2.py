import sys
import os
import re
import pandas as pd
from Bio import SeqIO
import argparse
import gzip


def parse_args():
    parser = argparse.ArgumentParser(description="Filter longest isoforms from a proteome.")
    parser.add_argument("input_proteome", help="Path to the input proteome file in faa.gz format")
    return parser.parse_args()




def determine_repo_and_mode(input_proteome):
    input_prefix = input_proteome.replace("protein.faa.gz", "")
    match = re.search(r'GCA_\d+|GCF_\d+', input_proteome)
    if not match:
        raise ValueError("Unknown repository type")

    identifier = match.group(0)
    if "GCA" in identifier:
        return "locus_tag", input_prefix
    elif "GCF" in identifier:
        return "GeneID", input_prefix




def filter_longest_isoforms(annotation, mode):
    proteins = annotation[(annotation['class'] == "with_protein") | (annotation['class'] == "protein_coding")]

    try:
        if re.search("-", proteins.loc[0,mode])[0]: # do not touch, it depends on the format of the annotation
            print("\nisoforms are annotated as second argument of", mode, "as in the example:", proteins.loc[0,mode],"\n")
            proteins[mode],proteins["isoform"]=proteins[mode].str.split("-",1).str
    except:
        print("regular annotation table format, having unique ids") 
        pass
    filtered = []
    genes = set(proteins[mode].dropna())
    for gene in genes:
        isoforms = proteins[proteins[mode] == gene]
        longest = isoforms[isoforms['product_length'] == isoforms['product_length'].max()]
        protein_id = longest['product_accession'].values[0]
        filtered.append(protein_id)

    return set(filtered)


def main():
    args = parse_args()

    mode, input_prefix = determine_repo_and_mode(args.input_proteome)
    print(input_prefix) 
    base_name = os.path.basename(input_prefix) 
    export_path = f"results/Proteomes/filtered_isoforms/{base_name}filtered.faa"
    print("... filtering", base_name)
    annotation = pd.read_csv(f"{input_prefix}feature_table.txt.gz", sep="\t")

    filtered_set = filter_longest_isoforms(annotation, mode)

    before_filt = 0
    after_filt = 0
    filtered_proteome = []

    with gzip.open(args.input_proteome, "rt") as proteome: # Note the "rt" mode
        for record in SeqIO.parse(proteome, 'fasta'):
            before_filt += 1
            if record.name in filtered_set:
                after_filt += 1
                filtered_proteome.append(record)

    SeqIO.write(filtered_proteome, export_path, 'fasta')

    print(f"Proteome total fasta records: {before_filt}")
    print(f"After the longest isoforms filtering: {after_filt}")


if __name__ == "__main__":
    main()
