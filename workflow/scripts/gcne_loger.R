library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

## Importing
args <- commandArgs(trailingOnly = TRUE)
Dmel_moultgene <- args[1]   # Workdir/dmel_moultgenes.csv  config
phylog =  args[2]    # Workdir/phylogeny_tmp.tsv
orthogroup = args[3]   #Results/loger/path2proteome_orthogroups.txt
output =  args[4]   # "Results/loger/loger_og_id_moult.tsv"


# genes of interest
dmel_moultgene_list <- readr::read_delim(Dmel_moultgene, delim = ",")
# phylogeny
phylogeny <- readr::read_delim(phylog, delim = "\t")
#not_timetree_orders <- read_delim("Workdir/excluded_orders.txt", delim="\t", col_names = F) # for cafe part, not needed now, need to come back to it later
# orthologer results files
loger_out <- read_table(orthogroup, col_names = FALSE, comment = "#") %>% as.data.frame()
names(loger_out) <- c("cluster_id", "ncbi_id", "cluster_type","seq_len", "seq_start", "seq_end","pid", "score","evalue")
loger_out <- separate(loger_out, 2, into = c("taxid", "ncbi_id"), sep="@", remove=TRUE)


## Filtering

# filtering out of orthologer results clusters of interest
# answer to the question: in which og do my genes of interest fall ?

loger_moult_full <- data.frame(matrix(nrow = 0, ncol = ncol(loger_out)))
cluster_ids <- c()
dmel_og <- filter(loger_out, taxid=="7227")

for (n in 1:nrow(dmel_moultgene_list)) {
  print(paste("searching the results for:", dmel_moultgene_list[n,"gene_name"]))
  reference <- dmel_moultgene_list[n,"NCBI_prot_ID"][[1,1]]
  target <- filter(dmel_og, ncbi_id == reference)
    if (nrow(target) != 0){
      print("FOUND")
      id <-target[1,"cluster_id"]
      if (id %in% cluster_ids == F){
        cluster_ids <- append(cluster_ids, id)
        og_df <- filter(loger_out, cluster_id == id)
        loger_moult_full <- rbind(loger_moult_full, og_df)
      }else if (id %in% cluster_ids){
        print(paste(reference, "falls in the previously filtered cluster", id))
        cluster_ids <- append(cluster_ids, id)}}} #implement
Dmel_moultgene <- mutate(dmel_moultgene_list, cluster_id = cluster_ids)

# retrieving genes specific for species/taxon other than dmel
others <- list(c("7070", "XP_970303.1", "cyp15a1", "sesquiterpenoid_pathway", "658858"))
#               c("7091", "NP_001140197.1", "cyp15c1", "sesquiterpenoid_pathway", "100286798"))

dig <- function(list_not_found){
  for (i in 1:length(list_not_found)) {
    animal <- list_not_found[[i]][1]
    gene <- list_not_found[[i]][2]
    print(paste("...retriving",gene, "in", animal))
    where <- filter(loger_out, taxid == animal, ncbi_id == gene)
    where <- select(where,cluster_id)[[1,1]]
        print(paste("the reference gene from this species falls in this cluster:", where))
        Dmel_moultgene <- add_row(Dmel_moultgene,
                                  Flybase_gene_ID = list_not_found[[i]][5],
                                  gene_name = list_not_found[[i]][3],
                                  known_function = list_not_found[[i]][4],
                                  NCBI_prot_ID = list_not_found[[i]][2],
                                  cluster_id = where)
        assign("Dmel_moultgene", Dmel_moultgene, envir = parent.frame())}}
dig(others)

# collapse rows from the same cluster
i <- 1
remove <- c()
while (i < nrow(Dmel_moultgene)){
  for (j in (i+1):nrow(Dmel_moultgene)){
    if (Dmel_moultgene[i,"cluster_id"] == Dmel_moultgene[j,"cluster_id"] && Dmel_moultgene[i,"known_function"] != Dmel_moultgene[j,"known_function"]) {
      remove <- c(remove, j)
      for (cols in 1:ncol(dmel_moultgene_list)){
        Dmel_moultgene[i,cols] <- base::paste(Dmel_moultgene[i,cols], Dmel_moultgene[j,cols], sep = ";")}}
    else if (Dmel_moultgene[i,"cluster_id"] == Dmel_moultgene[j,"cluster_id"] && Dmel_moultgene[i,"known_function"] == Dmel_moultgene[j,"known_function"]) {
      remove <- c(remove, j)
      for (cols in c(1:2,4:ncol(dmel_moultgene_list))){ Dmel_moultgene[i,cols] <- base::paste(Dmel_moultgene[i,cols], Dmel_moultgene[j,cols], sep = ";")}}}
  i <- i + 1} #implement
Dmel_moultgene <- Dmel_moultgene[-c(remove),]
rm(remove)


## Exporting result table having ogs from moulting genes of interest
result_table <- left_join(select(loger_moult_full, cluster_id, ncbi_id, taxid), select(Dmel_moultgene, gene_name, known_function, cluster_id), by = c("cluster_id" = "cluster_id"))

# rename for Valentine
new_result_table <- result_table %>% rename(protein_id = ncbi_id)
write.table(new_result_table, file = output, sep="\t", row.names = FALSE, quote = FALSE)


