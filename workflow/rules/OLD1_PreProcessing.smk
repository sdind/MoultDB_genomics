# The initial task refines the A3cat table by selecting high-quality genome assemblies, ensuring a balanced representation across species, and prepares data retrieval mechanisms for downstream protein analysis.



import os
import glob



rule a3cat_filtering:
# Parse A3cat table, filters assemblies that have annotation files, Selects the highest-quality assembly for each species, Excludes low-quality genomes and those with a limited number of protein-coding genes, downsamples overrepresented orders, Constructs FTP queries to retrieve protein.faa files for each species and finally gets the genome infos
    input:
        a3cat_table = config['a3cat_table']
    output:
        phylo = "results/a3cat/phylogeny_tmp.tsv",
        taxid = "results/a3cat/taxid4buscophile.tsv",
        genomes_map = "results/a3cat/genomes_map.tsv",
        features = "results/a3cat/feature_queries.csv",
        fasta_query = "results/a3cat/fasta_queries.csv"
    log:
        'logs/a3cat.log'
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
    Download proteome files and signal when the process is complete.
    '''
    input:
        features = rules.a3cat_filtering.output.features
    output:
        feat_out = "results/Proteomes/RawSeq/.featDone"
    log:
        'logs/get_feature.log'
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


rule get_proteome:
    '''
    Download proteome files, signal when the process is complete, and generate a list of proteome names.
    '''
    input:
        fasta_query = rules.a3cat_filtering.output.fasta_query
    output:
        Prot_out = "results/Proteomes/RawSeq/.ProtDone",
        proteome_list = "results/Proteomes/RawSeq/proteome_list.txt"
    log:
        'logs/get_proteome.log'
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
        touch {output.Prot_out}
        ls {params.outdir}/*_protein.faa.gz | sed 's/.*\\///;s/_protein.faa.gz//' > {output.proteome_list}
        """

rule filter_isoforms:
    '''
    Filter isoforms from the downloaded proteomes.
    '''
    input:
        rules.get_feature.output[0],
        proteome_list = "results/Proteomes/RawSeq/proteome_list.txt",
        proteome="results/Proteomes/RawSeq/{proteome_name}_protein.faa.gz"
    output:
        filtered="results/Proteomes/filtered_isoforms/{proteome_name}_filtered.faa"
    log:
        'logs/filter_isoforms/{proteome_name}.log'
    conda:
        '../envs/filter_isoforms.yaml'
    threads: 1
    resources:
        mem_mb = 10000
    params:
        runtime = '10:00:00',
        script = os.path.join(config['snakemake_dir_path'], 'workflow/scripts/filter_isoforms_v2.py')
    shell:
        """
        mkdir -p results/Proteomes/filtered_isoforms
        python {params.script} {input.proteome} > {output.filtered}
        """


rule all_filtered_proteomes:
    input:
        expand("results/Proteomes/filtered_isoforms/{proteome_name}_filtered.faa",
               proteome_name=glob_wildcards("results/Proteomes/RawSeq/{proteome_name}_protein.faa.gz").proteome_name)
    output:
        touch("results/Proteomes/filtered_isoforms/all_done.txt")


rule rename_proteomes:
    input:
        all_done="results/Proteomes/filtered_isoforms/all_done.txt",
        genomes_map="results/a3cat/genomes_map.tsv"
    output:
        touch("results/Proteomes/filtered_isoforms/rename_done.txt")
    log:
        "logs/filter_isoforms/rename_prot.log"
    shell:
        """
        bash -c '
        cd results/Proteomes/filtered_isoforms
        for genome in *faa ; do
            taxid=$(grep "$${genome:0:15}" {input.genomes_map} | cut -f 2)
            echo "... processing fasta labels: $$genome"
            sed "s/>/>$$taxid@/g" $$genome > $${taxid}.fa
            rm $$genome
        done
        touch {output}
        '
        """


rule create_orthologer_input:
    input:
        processed_files=rules.rename_proteomes.output[0]
    output:
        path2proteome="results/Proteomes/filtered_isoforms/path2proteome.tsv"
    log:
        "logs/interproscan/create_orthologer_input.log"
    shell:
        """
        cd results/Proteomes/filtered_isoforms
        path=$(pwd)
        lines=$(ls *fa | wc -l)
        ls -1 *fa > ext
        ls *fa | sed 's/\.[a-z]*$//g' > id
        printf "+\n%.0s" $(seq "$lines") > add
        printf "\t$path/\n%.0s" $(seq "$lines") > 2col
        paste add id -d "" > 1col
        paste 2col ext -d "" > 2col_full
        paste 1col 2col_full -d "" > "{output.path2proteome}"
        rm ext id add 1col 2col 2col_full
        """



rule run_interproscan_raw:
    input:
        proteome="results/Proteomes/RawSeq/{proteome_name}_protein.faa.gz"
    output:
        ipr_raw="results/Proteomes/interproscan_raw/{proteome_name}_ipr.tsv"
    conda:
        '../envs/InterProScan.yaml'
    threads: 8 
    resources:
        mem_mb = 150000
    log:
        "logs/interproscan/{proteome_name}.log"
    shell:
        """
        gunzip -c {input.proteome} > {input.proteome}.temp
        interproscan.sh -i {input.proteome}.temp -f TSV -o {output.ipr_raw} --cpu {threads}
        rm {input.proteome}.temp
        """


rule all_domains_proteomes:
    input:
        expand("results/Proteomes/interproscan_raw/{proteome_name}_ipr.tsv",
               proteome_name=glob_wildcards("results/Proteomes/RawSeq/{proteome_name}_protein.faa.gz").proteome_name)
    output:
        touch("results/Proteomes/interproscan_raw/domain_all_done.txt")




















# rule filter_isoforms:
#     '''
#     Filter isoforms from the downloaded proteomes.
#     '''
#     input:
#         os.path.join(os.path.dirname(rules.get_proteome_finished.output[0]), '{accession}.faa')
#     output:
#         filtered_fasta = "results/Proteomes/FilteredSeq/filt_{accession}.faa"
#     log:
#         'logs/filter_isoforms/filt_{accession}.log'
#     conda:
#         '../envs/filter_isoforms.yaml'
#     threads: 1
#     resources:
#         mem_mb = 5000
#     params:
#         runtime = '01:00:00'
#     shell:
#         "python workflow/scripts/filter_isoforms.py {input.fasta} > {output.filtered_fasta}"















# def get_fasta(wildcards):
#     checkpoints.get_proteome_finished.get()  # To make sure the DAG is re-evaluated to have all the result files
#     batch_input = rules.get_proteome_finished.output.outdir
#     batch_list = expand('results/Proteomes/RawSeq/{batch}.faa',
#                         batch=glob_wildcards(os.path.join(batch_input, '{batch}.faa')))
#     return batch_list




# rule filter_isoforms:
#     '''
#     Filter isoforms from the downloaded proteomes.
#     '''
#     input:
#         fasta = get_fasta
#     output:
#         filtered_fasta = directory("results/Proteomes/FilteredSeq/")
#     log:
#         'logs/filter_isoforms.log'
#     conda:
#         '../envs/filter_isoforms.yaml'
#     threads: 1
#     resources:
#         mem_mb = 5000
#     params:
#         runtime = '01:00:00'
#     shell:
#         "python workflow/scripts/filter_isoforms.py {input.fasta}"























# def get_sample_names(wildcards):
#     checkpoint_output = "results/Proteomes/RawSeq/"
#     filenames = glob.glob(checkpoint_output + "*.faa")
#     return [os.path.basename(f).rsplit('.', 1)[0] for f in filenames]





# import os
# import glob

# def get_sample_names(wildcards):
#     checkpoint_output = "results/Proteomes/RawSeq/"
#     filenames = glob.glob(checkpoint_output + "*.faa")
#     return [os.path.basename(f).rsplit('.', 1)[0] for f in filenames]




























# rule get_proteome: 
#     input:
#         features = rules.a3cat_filtering.output.features,
#         fasta_query = rules.a3cat_filtering.output.fasta_query
#     output:
#         outdir = directory("results/Proteomes/RawSeq/")
#     log:
#         'logs/get_proteome.log'
#     threads:1
#     resources:
#         mem_mb = 10000
#     params:
#         runtime = '20:00:00'
#     shell:
#         """
#         wget -i {input.features} -P {output.outdir} -w 3
#         wget -i {input.fasta_query} -P {output.outdir} -w 3
#         """















# rule aggregate_faa_files:
#     input:
#         directory("results/Proteomes/RawSeq/")
#     output:
#         "results/Proteomes/list_of_files.txt"
#     shell:
#         """
#         ls {input}/*.faa > {output}
#         """

# rule filter_isoforms:
#     input:
#         list_fa = "results/Proteomes/list_of_files.txt",
#         fa_dir = directory("results/Proteomes/RawSeq"),
#         fa = "results/Proteomes/RawSeq/{file}.faa"
#     output:
#         filtered_faa= "results/Proteomes/FilteredSeq/{file}_filtered.faa"
#     log:
#         "logs/filter_isoforms/{file}.log"
#     conda:
#         '../envs/filter_isoforms.yaml'
#     threads: 2
#     resources:
#         mem_mb = 5000
#     params:
#         runtime = '10:00:00'
#     shell:
#         """
#         python workflow/scripts/filter_isoforms.py {input.fa}  > {output.filtered_faa}"
#         """













# # checkpoint get_proteome_finished:
# #     input:
# #         directory("results/Proteomes/RawSeq/")
# #     output:
# #         directory("results/Proteomes/RawSeq/")
# #     shell:
# #         "touch {output}"


# # import os

# # def get_sample_names(wildcards):
# #     checkpoint_output = checkpoints.get_proteome_finished.get().output[0]
# #     filenames = glob.glob(checkpoint_output + "/*.faa")
# #     return [os.path.basename(f).rsplit('.', 1)[0] for f in filenames]



# # rule process_samples:
# #     input:
# #         get_sample_names
# #     output:
# #         expand("results/Proteomes/FilteredSeq/{sample}.filtered.faa", sample=get_sample_names())
# #     shell:
# #         "snakemake {output}"


# # rule filter_isoforms:
# #     input:
# #         fasta = "results/Proteomes/RawSeq/{sample}.faa"
# #     output:
# #         filtered_fasta = "results/Proteomes/FilteredSeq/{sample}.filtered.faa"
# #     log:
# #         'logs/filter_isoforms.{sample}.log'
# #     conda:
# #         '../envs/your_env.yaml'
# #     threads: 1
# #     resources:
# #         mem_mb = 5000
# #     params:
# #         runtime = '01:00:00'
# #     shell:
# #         "python workflow/scripts/filter_isoforms.py {input.fasta} > {output.filtered_faa}"








# rule aggregate_faa_files:
#     input:
#         directory("results/Proteomes/RawSeq/")
#     output:
#         "results/Proteomes/list_of_files.txt"
#     shell:
#         """
#         ls {input}/*.faa > {output}
#         """

# rule filter_isoforms:
#     input:
#         "results/Proteomes/list_of_files.txt",
#         "results/Proteomes/RawSeq/{file}.faa"
#     output:
#         "results/Proteomes/FilteredSeq/{file}_filtered.faa"
#     log:
#         "logs/filter_isoforms/{file}.log"
#     conda:
#         '../envs/filter_isoforms.yaml'
#     threads: 2
#     resources:
#         mem_mb = 5000
#     params:
#         runtime = '10:00:00'
#     shell:
#         """
#         python workflow/scripts/filter_isoforms.py {input[1]}
#         """

