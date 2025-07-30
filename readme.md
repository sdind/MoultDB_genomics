a3cat_filtered.tsv: Final genome information derived from the a3cat table, refined by selecting the highest-quality genome assemblies, downsampled overrepresented orders,
		    and ensured the inclusion of only those species for which proteomes were successfully downloaded, focusing on specific required columns.
                    Note to Valentine: Some columns were added and one column name was modified "Genbank Accession" to "RefSeq/GenBank Accession". This assembly number can be used to link genomes to all genes infos.


domains/"assembly_number"_"assembly_name"_filt_domains.tsv: Informations about protein domains of each genome/proteome. Retained only true positive InterPro domains.  
							    Note to Valentine: you can do the link between each file and its corresponding genome/protein thanks to the file name preset "assembly_number"

metadata_table/"assembly_number"_"assembly_name"_table.tsv: This dataset encompasses comprehensive metadata for all genes across all genomes/proteomes. 
							    The 'origin' column indicates whether the data was sourced from GenBank or RefSeq. For entries originating from
							    GenBank, there is an absence of both gene ID and gene name; in these cases, the 'locus_tag' is utilized instead.
							    Additionally, there is no 'transcript_id' for GenBank entries, necessitating the use of a URL for access.
							    Note to Valentine: Construct the URL using the following structure:: 
							     "https://www.ncbi.nlm.nih.gov/nuccore/"+"transcript_url_suffix"
 

prot_moult_pathways.tsv:  contains details on each identified moulting proteins, including their pathways, corresponding gene names, functions, and identifiers (= controlled 
			  voc), references. Giulia still need to complete the description 
			  Note to valentine: cluster_id correspond to orthougroups (probably not useful for you). 


orthogroups.tsv: contains orthogroups information about all the available proteins. 




