# Arthropod Moulting Gene Discovery Pipeline

This repository contains a Snakemake workflow for the selection, 
processing, and annotation of arthropod proteomes to identify orthologous 
groups associated with moulting pathways.

## Overview

The pipeline performs the following steps:

1. **Selection of High-Quality Proteomes**  
   Filters the [A3cat 
table](https://a3cat.unil.ch) to retain 
one representative genome per arthropod species, excluding low-quality 
assemblies and downsampling overrepresented lineages.

2. **Proteome Download and Isoform Filtering**  
   Downloads protein FASTA files and retains the longest isoform per 
gene.

3. **Metadata Extraction**  
   Extracts gene-level and protein-level metadata for downstream 
analysis.

4. **Orthologous Group Inference**  
   Runs [Orthologer](https://orthologer.ezlab.org) to infer 
orthologous groups across the selected species.

5. **Moulting Gene Identification**  
   Detects orthologs of known *Drosophila melanogaster* moulting genes, 
based on pathways curated by [Giulia Campli (PMID: 
39039636)](https://pubmed.ncbi.nlm.nih.gov/39039636/).

6. **Domain Annotation**  
   Annotates each filtered proteome using 
[InterProScan](https://www.ebi.ac.uk/interpro/), 
retaining high-confidence protein domains.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/)
- Conda (for environment management)
- Access to compute resources (recommended for InterProScan)

## Usage

Clone the repository:

```bash
git clone https://github.com/yourusername/MoultDB_genomics.git
cd MoultDB_genomics
```

Edit the config/config.yaml file to define input paths and parameters.

Run the pipeline:

```bash
snakemake --use-conda --cores 16
```

## Outputs

- a3cat_filtered.tsv: Final genome information derived from the a3cat table, refined by 
selecting the highest-quality genome assemblies, downsampled overrepresented orders, and ensured the inclusion of only those species for which proteomes were successfully downloaded, focusing on specific required columns.
                   
- domains/"assembly_number"_"assembly_name"_filt_domains.tsv: Informations about protein 
domains of each genome/proteome. Retained only true positive InterPro domains.  

- metadata_table/"assembly_number"_"assembly_name"_table.tsv: This dataset encompasses 
comprehensive metadata for all genes across all genomes/proteomes. The 'origin' column 
indicates whether the data was sourced from GenBank or RefSeq. For entries originating 
from GenBank, there is an absence of both gene ID and gene name; in these cases, the 
'locus_tag' is utilized instead. Additionally, there is no 'transcript_id' for GenBank entries, necessitating the use of a URL for access.

- prot_moult_pathways.tsv:  contains details on each identified moulting proteins, including their pathways, corresponding gene names, functions, and identifiers (= controlled voc), references. Giulia still need to complete the description 

- orthogroups.tsv: contains orthogroups information about all the available proteins. 




