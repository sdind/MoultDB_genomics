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


#### FOr now all I skip all this part, need to come back to it !! 

# ## Quantification of gene copy number per species (entire OG dataset)

# loger_counts <- count(loger_out, cluster_id, taxid)
# expand <- tidyr::expand(loger_counts,  cluster_id, phylogeny$TaxId)
# expand$key <- base::paste(expand$cluster_id, expand$`phylogeny$TaxId`, sep = "_") #considering species included in the analysis where putative orthologs were not detected (0 counts)
# loger_counts$key <- base::paste(loger_counts$cluster_id, loger_counts$taxid, sep = "_")
# loger_counts <- dplyr::left_join(expand, loger_counts, by = c("key" = "key"))%>%
#   dplyr::select(cluster_id.x, `phylogeny$TaxId`, n)%>%
#   tidyr::replace_na(list(n = 0L))
# loger_counts <- dplyr::rename(loger_counts, taxid = `phylogeny$TaxId`, cluster_id =cluster_id.x)
# loger_counts <- left_join(loger_counts, phylogeny, by = c("taxid" = "TaxId"))

# # median counts by order
# loger_median_counts <- loger_counts%>% group_by(cluster_id, Order) %>% summarise(median_count_by_order = median(n))
# large_variation <-  loger_median_counts %>% group_by(cluster_id) %>% summarise(max_min_diff = max(median_count_by_order)-min(median_count_by_order)) #cafe developers suggest to remove families with much variation between species

# ## Writing inputs for downstream analysis

# # cafe
# to_exclude <- filter(loger_median_counts,median_count_by_order >= 100)$cluster_id #cafe developers suggest to remove counts with at least 100 copies

# loger4cafe <- filter(loger_median_counts, !cluster_id %in% to_exclude)
# loger4cafe <- as.data.frame(tidyr::pivot_wider(loger4cafe, names_from = Order, values_from = median_count_by_order))
# loger4cafe <- select(loger4cafe, -not_timetree_orders$X1)
# loger4cafe <- mutate(loger4cafe, id =  paste("cluster",loger4cafe$cluster_id, sep="_")) %>%  relocate(id, .after = cluster_id)
# write.table(loger4cafe, file = "Workdir/loger4cafe.tsv", sep = "\t", row.names = FALSE, quote = F)

#aside_large_variation <- filter(large_variation, max_min_diff >= 60)$cluster_id
#loger4cafe_nolarge <- filter(loger4cafe, !cluster_id %in% aside_large_variation)
#write.table(loger4cafe_nolarge, file = "Workdir/loger4cafe_nolarge.tsv", sep = "\t", row.names = FALSE, quote = F)
