from Bio import SeqIO
import sys
import pandas as pd
import os 
import re

# the longest isoform is filtered : for each gene having the same locus_tag or GeneID (it depends on the repo of origin), the protein having the highest product length is taken

input_proteome=sys.argv[1]
input = input_proteome.replace("protein.faa.gz", "")

export="../filtered_isoforms/"+input+"filtered.faa"
if os.path.exists(export):
     print(export, "file already present in folder")
     quit()

input_feature=input+"feature_table.txt.gz"

print("... filtering", input_proteome)
repo=input[0:3] 
if repo == "GCA":
    mode="locus_tag" # do not touch, it depends on the format of the annotation
    print("Proteome has GenBank format")
elif repo == "GCF":
    mode="GeneID"
    print("Proteome has RefSeq format")

annotation=pd.read_csv(input_feature, sep="\t")
proteins=annotation[(annotation['class'] == "with_protein") | (annotation['class']  == "protein_coding")]

try:
    if re.search("-", proteins.loc[0,mode])[0]: # do not touch, it depends on the format of the annotation
        print("\nisoforms are annotated as second argument of", mode, "as in the example:", proteins.loc[0,mode],"\n")
        proteins[mode],proteins["isoform"]=proteins[mode].str.split("-",1).str
except:
    print("regular annotation table format, having unique ids") 
    pass

genes=set(list(proteins[mode]))
filtered=[]

for gene in genes:
    try:
        isoforms=proteins[(proteins[mode] == gene) & (proteins['class'] == "with_protein")]
        longest=isoforms[isoforms['product_length']==isoforms['product_length'].max()]
        protein_id=longest.loc[:,'product_accession'].values[0]
        filtered.append(protein_id)
    except:
        print("something went wrong, check feature table", mode, "and product_accession columns")
        pass

filtered_proteome=[]
before_filt=0
after_fil=0
with open(input_proteome, "r") as proteome: # implement
    records = SeqIO.parse(proteome, 'fasta')
    for record in records:
        name = record.name
        before_filt+=1
        if name in filtered:
            after_fil+=1
            filtered_proteome.append(record)

SeqIO.write(filtered_proteome, export, 'fasta')

print("proteome total fasta records:", before_filt, "\nafter the longest isoforms filtering:", after_fil, "\n")

#implement : don't do double for loop
