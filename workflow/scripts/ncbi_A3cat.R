#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
#library(HelpersMG)

## All available assemblies for Arthropoda are obtained from Arthropoda Assembly Assessment Catalogue (https://evofunvm.dcsr.unil.ch/table.html)
## and used to build strings to download data from NCBI FTP repositories

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]


## Importing catalogue data and cleaning
A3cat <- readr::read_delim(input_file, delim = "\t")
names(A3cat) <- base::gsub('\\s+', '_', names(A3cat))
A3cat <- mutate(A3cat, Annotation_Date = as.character(Annotation_Date))
A3cat$Assembly_Name <- base::gsub('\\s+', '_', A3cat$Assembly_Name)
A3cat$Species <- base::gsub('\\s+', '_', A3cat$Species)

## Data quick inspection 

# check how many unique species there are (despite n assemblies per species)
print(paste(nrow(A3cat), "assemblies"))
print(paste(length(unique(A3cat$Species)), "unique species"))

# filtering assemblies with annotation files
A3cat_anno <- dplyr::filter(A3cat, is.na(Annotation_Date) == F)
print(paste(nrow(A3cat_anno), "assemblies have annotation"))
print(paste(length(unique(A3cat_anno$Species)), "unique species"))

#choose if query only genomes with annotation (A3cat_anno) or anything (A3cat)
#data <- A3cat
data <- A3cat_anno

## Filtering best assemblies 

# filtering the best assembly for each species
data_dedup <- data.frame(matrix(ncol = ncol(data), nrow = 0))
names(data_dedup) <- names(data)
for (species in unique(data$Species)) { #using species instead of taxid prevents redundancy of strains (different organism but same taxid)
  print(paste("filtering duplicated assemblies with annotation for species:", species))
  find_assembly <- filter(data, Species == species) %>% arrange(desc(Arthropoda_Complete_BUSCO))
  find_assembly <- find_assembly[1,]
  data_dedup <- rbind(data_dedup, find_assembly)}

#filtering out low quality genome
thres <- 75 #choose threshold
data_dedup_qc <- dplyr::filter(data_dedup, Arthropoda_Complete_BUSCO >= thres) #NAs automatically excluded
print(paste(nrow(data_dedup_qc),"species have Arthropoda_Complete_BUSCO >", thres, ", percentage of non-insects is", round(nrow(filter(data_dedup_qc, Class != "Insecta"))/nrow(data_dedup_qc) * 100), "%"))

#filtering out assemblies with few protein-coding genes
print("filtering out assemblies with small number of protein coding genes:")
dplyr::filter(data_dedup_qc, `Protein-coding_genes` <= 8000) #to inspect what you are excluding
data_dedup_qc <- dplyr::filter(data_dedup_qc, `Protein-coding_genes` >= 8000)

## Downsizing the two most represented orders 

#downsampling most represented orders
print("Order distribution:")
order_distr <- data_dedup_qc %>% group_by(Order) %>% summarise(n_species=n()) %>% arrange(desc(n_species)) %>% as.data.frame() %>% print()
#write.table(order_distr, "order_distribution.tsv", sep = "\t", quote = F, row.names = F)

down <- function(df, rank){
  rank_only <- filter(df, Order == rank)
  down_rank <- data.frame(matrix(ncol= ncol(rank_only), nrow = 0))
  names(down_rank) <- names(rank_only)
  for (family in unique(rank_only$Family)){
  if ( family %in% c("Tenebrionidae","Drosophilidae") == F){
    best <- filter(rank_only, Family == family) %>% as.data.frame()
    set <-paste(rank, "_Complete_BUSCO", sep="")
    if (set %in% names(best) == F) { set <- "Arthropoda_Complete_BUSCO"}
    best <- best[order(best[,set], decreasing = TRUE),]
    best <- best[1,] #taking only one best assembly for each FAMILY
    #if (nrow(best) > 1){best <- best[1:2,]} else {best <- best[1,]}  #taking 2 assembly for each family
    down_rank <- rbind(down_rank, best)
    if (nrow(down_rank) >= 15){
      print(paste(rank, "has reached 15 species"))
      down_rank <- down_rank[order(down_rank[,set], decreasing = TRUE),]
      down_rank <- down_rank[1:15,]}
 } }
    return(assign(paste("down_", rank, sep=""), down_rank, envir = parent.frame()))
  #return(print(down_rank))
}

resized <- data.frame(matrix(ncol=ncol(data_dedup_qc), nrow = 0))
for (insect in unique(select(filter(data_dedup_qc, Class == "Insecta"), Order))$Order){
 if( filter(order_distr, Order == insect)$n_species > 20){
  print(paste("downsampling", insect))
  resized <- rbind(resized, down(data_dedup_qc, insect))
 }
 }
                             
genomes2query <- rbind(resized, filter(data_dedup_qc, Assembly_Name %in% c("Release_6_plus_ISO1_MT","Tcas5.2")), filter(data_dedup_qc, !Order %in%  unique(resized$Order)))
genomes2query %>% group_by(Order) %>% summarise(n_species=n()) %>% arrange(desc(n_species)) %>% as.data.frame() %>% print()                             

## Building ftp queries for protein.faa files for each species

print("building ncbi ftp queries for proteome fasta files")

queries <- c()
table <- c()
assemblies <- c()
refseq_counter <- 0
assembly_taxid_phylo <- data.frame(matrix(nrow = 0, ncol = 7))

for (i in 1:nrow(genomes2query)) {
  assembly_name <- genomes2query[i,"Assembly_Name"]
  species_name <- genomes2query[i,"Species"]
  if (is.na(genomes2query[i,"Refseq_Accession"])) {
    assembly_number <- genomes2query[i,"Genbank_Accession"] #takes Genbank if no RefSeq repo is available
    assembly_taxid_phylo  <- filter(genomes2query, Genbank_Accession == assembly_number) %>% 
      dplyr::select(TaxId, Species, SubPhylum, Class, Order, Family, Genus) %>%data.frame(accession = assembly_number) %>% rbind(assembly_taxid_phylo) 
    queries <- base::append(queries, paste("ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/", species_name, "/all_assembly_versions/", assembly_number,"_", assembly_name, "/", assembly_number,"_", assembly_name, "_protein.faa.gz", sep = ""))
    table <- base::append(table, paste("ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/", species_name, "/all_assembly_versions/", assembly_number,"_", assembly_name, "/", assembly_number,"_", assembly_name, "_feature_table.txt.gz", sep = ""))
  } else {
    assembly_number <- genomes2query[i, "Refseq_Accession"] 
  assembly_taxid_phylo  <- filter(genomes2query, Refseq_Accession == assembly_number) %>% 
    dplyr::select(TaxId, Species, SubPhylum, Class, Order, Family, Genus) %>% data.frame(accession = assembly_number)%>% rbind(assembly_taxid_phylo)
  queries <- base::append(queries, paste("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/", species_name, "/all_assembly_versions/", assembly_number,"_", assembly_name, "/", assembly_number,"_", assembly_name, "_protein.faa.gz", sep = ""))
  table <- base::append(table, paste("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/", species_name, "/all_assembly_versions/", assembly_number,"_", assembly_name, "/", assembly_number,"_", assembly_name, "_feature_table.txt.gz", sep = ""))
  #refseq_counter <- refseq_counter +1
  }
}

#add here path to species of interest already known to be missed
print("adding oncopeltus fasciatus")
assembly_taxid_phylo <- tibble::add_row(assembly_taxid_phylo, 
                                        TaxId = 7536, Species = "Oncopeltus_fasciatus", 
                                        Class = "Insecta", Order = "Hemiptera", Family = "Lygaeidae", Genus = "Oncopeltus", SubPhylum = "Hexapoda",
                                        accession = "GCA_000696205.1")
queries <- append(queries, "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/Oncopeltus_fasciatus/all_assembly_versions/GCA_000696205.1_Ofas_1.0/GCA_000696205.1_Ofas_1.0_protein.faa.gz")
table <- append(table, "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/Oncopeltus_fasciatus/all_assembly_versions/GCA_000696205.1_Ofas_1.0/GCA_000696205.1_Ofas_1.0_feature_table.txt.gz")

#check
head(queries)
refseq_counter <-  which(is.na(genomes2query$Refseq_Accession) == F) %>% length()
print(paste("total queries:", length(queries)))
print(paste("percentage of proteomes fetched from RefSeq:", round((refseq_counter/length(queries))*100), "%"))


## Writing temporary files for downstream analysis

print("writing temporary files")
write.table(assembly_taxid_phylo, file = "results/a3cat/phylogeny_tmp.tsv", row.names = F, col.names = T, quote = F, sep = "\t") #as finder_phylogeny.tsv
write.table(dplyr::select(assembly_taxid_phylo, TaxId), file = "results/a3cat/taxid4buscophile.tsv", row.names = F, col.names = F, quote = F, sep = "\t") #input for buscophile
write.table(dplyr::select(assembly_taxid_phylo, accession, TaxId, Species), file = "results/a3cat/genomes_map.tsv", row.names = F, col.names = F, quote = F, sep = "\t" ) # mappings genome-taxid-species for fasta relabelling

## querying

#getting genomes
#print("Querying ncbi ftps ... ")
utils::write.table(queries, file ="results/a3cat/fasta_queries.csv", row.names = F, col.names = F, quote = F, sep = ",") 
utils::write.table(table, file ="results/a3cat/feature_queries.csv", row.names = F, col.names = F, quote = F, sep = ",") 
#HelpersMG::wget(queries)
