library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

## Importing
args <- commandArgs(trailingOnly = TRUE)
Dmel_moultgene <- args[1]   # Workdir/Dmel_moultgene_allIDs_R_narrow_list3.csv
phylog =  args[2]    # Workdir/phylogeny_tmp.tsv
orthogroup = args[3]   #Results/loger/path2proteome_orthogroups.txt
output =  args[4]   # "Results/loger/loger_og_id_moult.tsv"


# genes of interest
Dmel_moultgene_allIDs_narrow_list <- readr::read_delim(Dmel_moultgene, delim = ",") %>% as.data.frame()
# phylogeny
phylogeny <- readr::read_delim(phylog, delim = "\t")
# orthologer results files
loger_out <- read_table2(orthogroup, col_names = FALSE, comment = "#") %>% as.data.frame()#path and name will be adjusted
names(loger_out) <- c("cluster_id", "ncbi_id", "cluster_type","seq_len", "seq_start", "seq_end","pid", "score","evalue")

loger_out <- mutate(loger_out, taxid = str_extract(ncbi_id, "[:digit:]{3,}(?=@)")) #preprocessed fasta have taxid for each seq

## Filtering

# filtering out of orthologer results clusters of interest
# answer to the question: in which og do my genes of interest fall ?

loger_moult_full <- data.frame(matrix(nrow = 0, ncol = ncol(loger_out)))
cluster_ids <- c()

dmel_og <- filter(loger_out, taxid=="7227")

for (n in 1:nrow(Dmel_moultgene_allIDs_narrow_list)) {
  print(paste("searching the results for:", Dmel_moultgene_allIDs_narrow_list[n,"gene_name"]))
  A <- as.list(base::unlist(base::strsplit(Dmel_moultgene_allIDs_narrow_list[[n,"NCBI_prot_ID"]], ";")))
  for (m in 1:nrow(dmel_og)) {
    B <- str_remove(dmel_og[m,"ncbi_id"],"7227@" )
    if (length(intersect(A, B)) != 0){
      print("FOUND")
      id <- dmel_og[m,"cluster_id"]
      if (id %in% cluster_ids == F){
        cluster_ids <- append(cluster_ids, id)
        og_df <- filter(loger_out, cluster_id == id)
        loger_moult_full <- rbind(loger_moult_full, og_df)
      }else if (id %in% cluster_ids){
        print(paste(B, "falls in the previously filtered cluster", id))
        cluster_ids <- append(cluster_ids, id)
      }
    }
  }
} #implement

Dmel_moultgene <- mutate(Dmel_moultgene_allIDs_narrow_list, cluster_id = cluster_ids)

# collapse rows from the same cluster
i <- 1
remove <- c()
while (i < nrow(Dmel_moultgene)){
  for (j in (i+1):nrow(Dmel_moultgene)){
    if (Dmel_moultgene[i,"cluster_id"] == Dmel_moultgene[j,"cluster_id"] && Dmel_moultgene[i,"known_function"] != Dmel_moultgene[j,"known_function"]) {
      remove <- c(remove, j)
      for (cols in 1:ncol(Dmel_moultgene_allIDs_narrow_list)){ Dmel_moultgene[i,cols] <- base::paste(Dmel_moultgene[i,cols], Dmel_moultgene[j,cols], sep = ";")}
    }
    else if (Dmel_moultgene[i,"cluster_id"] == Dmel_moultgene[j,"cluster_id"] && Dmel_moultgene[i,"known_function"] == Dmel_moultgene[j,"known_function"]) {
      remove <- c(remove, j)
      for (cols in c(1:4,6:ncol(Dmel_moultgene_allIDs_narrow_list))){ Dmel_moultgene[i,cols] <- base::paste(Dmel_moultgene[i,cols], Dmel_moultgene[j,cols], sep = ";")
      }
    }
  }
  i <- i + 1
} #implement

Dmel_moultgene <- Dmel_moultgene[-c(remove),]
rm(remove)

loger_moult <- loger_moult_full %>%
  select(cluster_id, ncbi_id, cluster_type, taxid)

## Exporting result table having ogs from moulting genes of interest
write.table(left_join(select(loger_moult, cluster_id, ncbi_id, taxid), select(Dmel_moultgene, gene_name, known_function, cluster_id), by = c("cluster_id" = "cluster_id")), file = output, sep="\t", row.names = FALSE, quote = F)

