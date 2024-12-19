library(tidyverse)
library(DECIPHER)
library(here)
library(GrpString)
library(purrr)
library(optparse)

#Arguments
option_list = list(
  make_option(c("-t", "--path_tax"), action="store", default=NA, type='character',
              help="path to taxonomy file in rCRUX format"),
  
  make_option(c("-f", "--path_fasta"), action="store", default=NA, type='character',
              help="path to fasta file in rCRUX format"),
  
  make_option(c("-b", "--barcode_name"), action="store", default=NA, type='character',
              help="name of target marker gene or barcode"), 
  make_option(c("-o", "--output_dir"), action="store", default=NA, type='character',
              help="path to output directory"),
  make_option(c("-l", "--barcode_length"), action="store", default=NA, type='numeric',
              help="average length of marker gene, e.g. 186 [MiFish 12S =]")
)


opt = parse_args(OptionParser(option_list=option_list))

tax <- read.table(opt$path_tax, sep="\t", header=F, fill=FALSE) %>% dplyr::rename(Accession = V1,taxonomic_path=V2) 
fasta <- readDNAStringSet(opt$path_fasta) 

## All Fishes

tax %>% 
  unite(., `Accession:taxonomic_path`, c("Accession","taxonomic_path"), sep = ":", remove=F) -> tax_e

as.data.frame(fasta) %>% as_tibble(rownames = "Accession") %>% 
  left_join(tax_e) -> fasta_df

fasta_df %>% 
  filter(., !is.na(taxonomic_path)) %>% 
  filter(., !str_detect(taxonomic_path,"NA;NA;NA;NA")) %>% 
  filter(., !str_detect(taxonomic_path,"Bacteria;")) %>%
  filter(., !str_detect(taxonomic_path,"Archaea;")) -> fasta_df_clean

fasta[fasta_df_clean$Accession] -> fasta_clean

names(fasta_clean) <- fasta_df_clean$`Accession:taxonomic_path`

cutoff_score = 1/opt$barcode_length

clusters <- DECIPHER::Clusterize(fasta_clean, cutoff=cutoff_score, invertCenters=FALSE)

apply(clusters, 2, function(x) max(abs(x))) # number of clusters

clusters %>% 
  rownames_to_column(., var = "Name") %>% 
  arrange(cluster) %>% 
  separate(Name, into =c("Accession","Name"), sep = ":", extra = "merge") %>% 
  mutate(., Name=str_remove_all(Name, "\\_[0-9]+") ) %>% 
  dplyr::select(cluster,Name,Accession) -> cluster_IDS

cluster_IDS %>% 
  dplyr::select(-Accession) %>% 
  distinct()  %>% 
  filter(., !is.na(Name)) %>% 
  group_by(cluster) %>% summarise(n=sum(n())) %>% filter(n>1) -> duplicate_reference_cluster

cluster_IDS %>% 
  filter(., !is.na(Name)) %>% 
  filter(., !cluster %in% duplicate_reference_cluster$cluster) %>% 
  dplyr::select(cluster, lca_taxonomic_path = Name) %>% 
  mutate(., data=NA) %>% 
  group_by(cluster,lca_taxonomic_path) %>% 
  nest()-> trust_worthy_IDS

cluster_IDS %>% 
  group_by(cluster) %>% 
  nest(Name,Accession) %>% 
  dplyr::rename(Accesions = data) -> accession_info_cluster_IDS

trust_worthy_IDS %>% 
  left_join(accession_info_cluster_IDS) -> trust_worthy_IDS_accessions

cluster_IDS %>% 
  filter(., cluster %in% duplicate_reference_cluster$cluster)-> sequence_clusters_multiple_names

sequence_clusters_multiple_names %>% 
  dplyr::rename(taxonomic_path=Name) %>% 
  dplyr::select(-Accession) %>% 
  distinct() %>% 
  group_by(cluster) %>% 
  nest() %>% 
  mutate(data = purrr::map(data,  function(.x){
    .x %>%
      ungroup() %>% 
      mutate(idx= match(taxonomic_path, unique(taxonomic_path))) %>% 
      separate(taxonomic_path, into=c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";", remove=F) %>% 
      pivot_wider(names_from = idx, values_from = c(taxonomic_path, Domain,Phylum,Order,Class, Family, Genus, Species),  names_glue = "{.value}:{idx}" ) %>% 
      rowwise %>% 
      mutate(Same_Species = n_distinct(unlist(across(starts_with('Species'), 
                                                     ~ as.character(.x)))) == 1,
             Same_Genus = n_distinct(unlist(across(starts_with('Genus'), 
                                                   ~ as.character(.x)))) == 1,
             Same_Family = n_distinct(unlist(across(starts_with('Family'), 
                                                    ~ as.character(.x)))) == 1,
             Same_Class = n_distinct(unlist(across(starts_with('Order'), 
                                                   ~ as.character(.x)))) == 1,
             Same_Order = n_distinct(unlist(across(starts_with('Class'), 
                                                   ~ as.character(.x)))) == 1,
             Same_Phylum = n_distinct(unlist(across(starts_with('Phylum'), 
                                                    ~ as.character(.x)))) == 1,
             Same_Domain = n_distinct(unlist(across(starts_with('Domain'), 
                                                    ~ as.character(.x)))) == 1) %>%
      ungroup() %>% 
      mutate(., lca_taxonomic_path = case_when(Same_Species==TRUE ~ `taxonomic_path:1`,
                                               Same_Domain==FALSE ~ paste0(";",";",";",";",";",";"),
                                               Same_Phylum==FALSE ~ paste0(`Domain:1`,";",";",";",";",";",";"),
                                               Same_Class==FALSE ~ paste0(`Domain:1`,";",`Phylum:1`,";",";",";",";",";"),
                                               Same_Order==FALSE ~ paste0(`Domain:1`,";",`Phylum:1`,";",`Class:1`,";",";",";",";"),
                                               Same_Family==FALSE ~ paste0(`Domain:1`,";",`Phylum:1`,";",`Class:1`,";",`Order:1`,";",";",";"),
                                               Same_Genus==FALSE ~ paste0(`Domain:1`,";",`Phylum:1`,";",`Class:1`,";",`Order:1`,";",`Family:1`,";",";"),
                                               Same_Species==FALSE ~ paste0(`Domain:1`,";",`Phylum:1`,";",`Class:1`,";",`Order:1`,";",`Family:1`,";",`Genus:1`,";"),
                                               TRUE ~"weird")) })) %>% 
  unnest() %>% 
  group_by(cluster,lca_taxonomic_path) %>% 
  nest() -> sequence_clusters_multiple_names_nested

sequence_clusters_multiple_names %>% 
  group_by(cluster) %>% 
  nest(Name,Accession) %>% 
  dplyr::rename(Accesions = data) -> accession_info_multiple_names

sequence_clusters_multiple_names_nested %>% 
  left_join(accession_info_multiple_names) -> sequence_clusters_multiple_names_nested_accessions


rbind(trust_worthy_IDS_accessions,sequence_clusters_multiple_names_nested_accessions) %>% 
  arrange(cluster)-> full_merged_data_set

saveRDS(full_merged_data_set, here(paste0(opt$output_dir,opt$barcode_name,"_full_merged_data_set.RDS")))

full_merged_data_set %>%
  filter(., str_detect(lca_taxonomic_path, "Eukaryota;Chordata;" )) %>% 
  separate(lca_taxonomic_path, into=c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";", remove=F) -> fish_separated

saveRDS(fish_separated, here(paste0(opt$output_dir,opt$barcode_name,"_vertebrates_only_merged_data_set.RDS")))


