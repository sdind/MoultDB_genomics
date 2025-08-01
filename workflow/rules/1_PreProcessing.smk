import os
import glob


rule a3cat_filtering:
    """
    This rule refines the A3cat table by selecting the highest-quality genome assemblies.
    It parses the table, selects the highest-quality assembly per species, excludes low-quality genomes with limited protein-coding genes, downsamples overrepresented orders, and constructs FTP queries to retrieve protein fasta files and to get their genome infos
    """
    input:
        a3cat_table = config['a3cat_table']
    output:
        phylo = "results/a3cat/phylogeny_tmp.tsv",
        taxid = "results/a3cat/taxid4buscophile.tsv",
        genomes_map = "results/a3cat/genomes_map.tsv",
        features = "results/a3cat/feature_queries.csv",
        fasta_query = "results/a3cat/fasta_queries.csv"
    log:
        'logs/a3cat_filt/a3cat.log'
    conda:
        '../envs/a3cat.yaml'
    threads:1
    resources:
        mem_mb = 5000
    params:
        runtime = '01:00:00'
    shell:
        "Rscript workflow/scripts/ncbi_A3cat.R {input.a3cat_table}"




rule get_feature:
    '''
    Downloads proteome files using the feature queries produced by the a3cat_filtering rule. This rule reads a list of FTP links to download and store proteome files in a designated directory.
    '''
    input:
        features = rules.a3cat_filtering.output.features
    output:
        feat_out = "results/Proteomes/RawSeq/.featDone"
    log:
        'logs/get_feature/get_feature.log'
    threads: 1
    resources:
        mem_mb = 10000
    params:
        runtime = '10:00:00',
        outdir = directory("results/Proteomes/RawSeq")

    shell:
        """
        while read line; do
            wget $line -P {params.outdir} -w 20 || echo "Failed to download $line" >> {log}
        done < {input.features}
        touch {output.feat_out}
        """



# Checkpoint to collect proteome names in order to execute later rules (as we don't know the proteome names a priori!)
checkpoint get_proteome:
    '''
    A checkpoint rule that downloads proteome fasta files based on queries from the a3cat_filtering rule. This rule also generates a list of proteome names for downstream processing.
    '''
    input:
        fasta_query = "results/a3cat/fasta_queries.csv"
    output:
        proteome_list = "results/Proteomes/RawSeq/proteome_list.txt"
    log:
        "logs/get_proteome/get_proteome.log"
    threads: 2
    resources:
        mem_mb = 30000
    params:
        runtime = '10:00:00',
        outdir = directory("results/Proteomes/RawSeq")
    shell:
        """
        mkdir -p {params.outdir}
        while read line; do
            wget $line -P {params.outdir} -w 20 || echo "Failed to download $line" >> {log}
        done < {input.fasta_query}
        ls {params.outdir}/*_protein.faa.gz | sed 's/.*\///;s/_protein.faa.gz//' > {output.proteome_list}
        """




rule filter_table:
   '''
    Refines the a3cat table post-proteome download, ensuring inclusion of only those species with successfully downloaded proteomes and only specific required columns. 
    '''
    input:
        a3cat_table = config['a3cat_table'],
        proteome_list = "results/Proteomes/RawSeq/proteome_list.txt"
    output:
        selected_data = "results/moultDB_input/a3cat_filtered.tsv"
    conda:
        '../envs/filter_isoforms.yaml'
    log:
        "logs/moultDB_input/filter_and_select.log"
    script:
        "../scripts/filter_a3cat.py"





rule filter_isoforms:
    '''
    This rule processes each proteome file to retain only the longest isoform for each gene. 
    '''
    input:
        rules.get_feature.output[0],
        proteome="results/Proteomes/RawSeq/{proteome_name}_protein.faa.gz"
    output:
        filtered="results/Proteomes/filtered_isoforms/{proteome_name}_filtered.faa"
    log:
        "logs/filter_isoforms/filter_isoforms_{proteome_name}.log"
    conda:
        '../envs/filter_isoforms.yaml'
    threads: 4
    resources:
        mem_mb = 10000
    params:
        runtime = '35:00:00',
        script="workflow/scripts/filter_isoforms_v2.py"
    shell:
        """
        mkdir -p results/Proteomes/filtered_isoforms
        python {params.script} {input.proteome}
        """




rule retrieve_metadata:
    """
    Extracts and compiles metadata for each protein within the filtered proteomes and create a table with all these info. 
    Metadata includes gene ID, gene description, protein_length, CDS_protein_Id, CDS_protein_Id-Xref (for url). 
    """
    input:
        filtered = "results/Proteomes/filtered_isoforms/{proteome_name}_filtered.faa"
    output:
        metadata_table = "results/Proteomes/metadata_table/{proteome_name}_table.tsv"
    log:
        "logs/Proteomes/metadata_table/metadata_{proteome_name}.log"
    threads: 2
    conda:
        '../envs/retrieve_metadata.yaml'
    params:
        runtime = '70:00:00'
    resources:
        mem_mb = 10000
    script:
        "../scripts/retrieve_prot_metadata.py"




def get_filtered_proteome_names(wildcards):
    checkpoint_output = checkpoints.get_proteome.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        proteome_names = [line.strip() for line in f]
    return expand("results/Proteomes/filtered_isoforms/{proteome_name}_filtered.faa", 
                  proteome_name=proteome_names)




rule process_fasta_labels:
    '''
    Renames fasta file names and seqeunces headers with taxonomy IDs instead of proteome names to meet orthologer requirements.
    '''
    input:
        filtered_proteomes = get_filtered_proteome_names,   # just to be sure that they all finished before processing 
        genomes_map = "results/a3cat/genomes_map.tsv"
    output:
        touch("results/Proteomes/filtered_isoforms/processed_fasta_label.done")
    log:
        "logs/process_fasta_labels.log"
    threads: 1
    resources:
        mem_mb = 1000
    params:
        runtime = '03:00:00',
        filt_dir = "results/Proteomes/filtered_isoforms"
    shell:
        """
        cd {params.filt_dir}
        for genome in *faa; do
            taxid=$(grep ${{genome:0:15}} ../../../{input.genomes_map} | cut -f 2)
            sed "s/>/>${{taxid}}@/g" ${{genome}} > ${{taxid}}.fa
            #rm ${{genome}}
        done
        """
## to change later: genome:0:15 not robust (Giulia's command), to change with regex + delete all files when workflow is done



rule generate_orthologer_input:
    '''
    Prepares input for Orthologer by creating a path2proteome.tsv file, which lists the paths to fasta files in the format required by Orthologer (need to start with a +)
    '''
    input:
        processed_fa_labels = "results/Proteomes/filtered_isoforms/processed_fasta_label.done"
    output:
        orthologer_input = "results/a3cat/path2proteome.tsv"
    log:
        "logs/path2proteome/generate_orthologer_input.log"
    params:
        runtime = '00:05:00',
        filt_dir = "results/Proteomes/filtered_isoforms"
    shell:
        """
        cd {params.filt_dir}
        path=$(pwd)
        > ../../../{output.orthologer_input}   #Empty the output file if it already exists
        for file in *fa; do
            filename=$(basename "$file" .fa)
            echo "+${{filename}}\t${{path}}/${{file}}" >> ../../../{output.orthologer_input}
        done
        """




rule run_orthologer_S1:
    """
    Executes the Orthologer to identify orthologous groups across multiple species. 
    This rule automates Orthologer setup and execution using run_loger_prompt.sh to handle interactive prompts.
    """
    input:
        path2proteome = "results/a3cat/path2proteome.tsv"
    output:
        OGs = "results/orthologer/mydata/Results/path2proteome_orthogroups.txt"
    log:
        "logs/orthologer/orthologerS1_run.log"
    conda:
        '../envs/orthologer_040923.yaml'
    params:
        orthologer_dir = "/work/FAC/FBM/DEE/mrobinso/moult/giulia/orthologer_v3/orthologer_3.0.2/ORTHOLOGER-3.0.2/bin",
        script_prompt = os.path.join(config['snakemake_dir_path'],'workflow/scripts/run_loger_prompt.sh'),
        runtime = '10:10:00'
    threads: 2
    resources:
        mem_mb = 500
    shell:
        """
        mkdir -p results/orthologer/mydata
        cd results/orthologer/mydata
        cp ../../../{input.path2proteome} .
        {params.orthologer_dir}/manage_project.sh -s -f path2proteome.tsv
        {params.script_prompt}
        ./setup_project_sdind2.sh     # to change in the future to be more robust !!
        cp /work/FAC/FBM/DEE/mrobinso/moult/giulia/orthologer_v3/orthologer_3.0.2/ORTHOLOGER-3.0.2/orthologer_conf.sh .
        PYTHON3_PATH=$(which python3)
        sed -i "s|PYTHON3=\"/work/FAC/FBM/DEE/mrobinso/moult/giulia/mambaforge/envs/orthologer_040923/bin/python3\"|PYTHON3=\"$PYTHON3_PATH\"|" orthologer_conf.sh 
        ./orthologer.sh -xp
        ./orthologer.sh -r all
        """


rule reformat_OG_out:
    input:
        OGs = "results/orthologer/mydata/Results/path2proteome_orthogroups.txt"
    output:
        OGs_reformat = "results/moultDB_input/orthogroups.tsv"
    shell:
        r"""
        echo -e "orthogroup_id\ttaxid\tprotein_id\tversion" > {output.OGs_reformat};
        awk -v date="$(date '+%Y-%m-%d')" 'BEGIN{{OFS="\t"}} !/^#/{{split($2, a, "@"); print $1, a[1], a[2], date}}' {input.OGs} >> {output.OGs_reformat}
        """



rule run_filter_OGs:
    """
    Filters the orthogroups to retain only moulting genes in order to attribute the corresponding moulting pathways described by Giulia Camplu. 
    This rule uses a predefined carefully curated list of Drosophila melanogaster moulting gene IDs as reference 
    """
    input:
        Dmel_moultgene = config["Dmel_moultgeneIDs"],
        pyhlo_tsv = config["phylogeny_tsv"],
        OGs = "results/orthologer/mydata/Results/path2proteome_orthogroups.txt" 
    output:
        filtered_OGs_pathways = "results/orthologer/mydata/Results/loger_og_id_moult.tsv"
    log:
        "logs/orthologer/filter_OGs.log"
    conda:
        '../envs/filter_OGs.yaml'
    params:
        script_gcne = "workflow/scripts/gcne_loger.R",
        runtime = '10:10:00'
    threads: 2
    resources:
        mem_mb = 10000
    shell:
        "Rscript {params.script_gcne} {input.Dmel_moultgene} {input.pyhlo_tsv} {input.OGs} {output.filtered_OGs_pathways}"



rule add_genome_ids_and_pathway_info:
    input:
        filtered_OGs_pathways = "results/orthologer/mydata/Results/loger_og_id_moult.tsv",
        genomes_map = "results/a3cat/genomes_map.tsv",
        path_voc = config["pathway_controlled_voc"]
    output:
        pathways_output = "results/moultDB_input/prot_moult_pathways.tsv"
    log:
        "logs/moultDB_input/add_genome_ids_and_pathway_info.log"
    params:
        runtime = '00:10:00'
    threads: 1
    resources:
        mem_mb = 2000
    shell:
        """
        # Convert path_voc.tsv to a tab-separated file for easier processing
        awk 'BEGIN{{FS=","; OFS="\t"}} {{print $1,$2,$3,$4}}' {input.path_voc} > path_voc_tab.tsv
        
        # Map genome IDs and append pathway information based on known_function match
        awk 'BEGIN{{FS=OFS="\t"}} 
            NR==FNR {{genome_map[$2]=$1; next}} 
            FNR==1 && NR!=FNR {{print $0, "genome_id", "identifier", "reference", "description"; next}}
            NR!=FNR {{pathway_info[$1]=$2 FS $3 FS $4}} 
            NR==FNR || FNR>1 {{print $0, genome_map[$3], pathway_info[$5]}}' {input.genomes_map} path_voc_tab.tsv {input.filtered_OGs_pathways} > {output.pathways_output}
        """





rule run_interproscan_raw:
    '''
    Executes InterProScan on each filtered proteome to identify protein domains, leveraging parallel processing for efficiency. 
    '''
    input:
        filtered = "results/Proteomes/filtered_isoforms/{proteome_name}_filtered.faa"
    output:
        domains = "results/Proteomes/interproscan/domains/{proteome_name}.domains"
    log:
        "logs/interproscan/{proteome_name}_domains.log"
    params:
        output_dir=directory("results/Proteomes/interproscan/domains"),
        runtime = '70:00:00'
    threads: 8
    resources:
        mem_mb = 95000
    singularity:
        'docker://interpro/interproscan:5.65-97.0'
    shell:
        """
        mkdir -p {params.output_dir}
        /opt/interproscan/interproscan.sh -dp -i {input.filtered} -f TSV -o {output.domains} --cpu {threads}
        """



rule filt_domains:
    '''
    Filters the domain results to retain only true positive InterPro domains. Store IPR_id, start, end in a output table.
    '''
    input:
        domains = "results/Proteomes/interproscan/domains/{proteome_name}.domains"
    output:
        filtered_domains = "results/moultDB_input/domains/{proteome_name}_filt_domains.tsv"
    log:
        "logs/Proteomes/filt_domains/{proteome_name}_filt_domains.log"
    threads: 2
    params:
        runtime = '25:00:00'
    resources:
        mem_mb = 120000
    script:
        "../scripts/filt_domains.py"















