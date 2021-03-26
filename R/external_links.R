#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
generate_external_links <- function(names, activities){
  
}



.cmpd_links <- function(names){
  db.names <- readRDS(syn$get("syn12978912")$path)
  
  chembl.internal.ids <- db.names %>% 
    filter(database == "chembl") %>%
    select(external_id, internal_id, database) %>% 
    distinct()
  
  chembl.links <- read.table(syn$get("syn12972665")$path, sep = "\t", quote = "", comment.char = "", header = T) %>% 
    select(molregno, chembl_id) %>% 
    distinct() %>% 
    set_names(c("external_id", "link")) %>% 
    mutate(external_id = as.character(external_id)) %>% 
    right_join(chembl.internal.ids) %>% 
    mutate(link = paste0("<a href='https://www.ebi.ac.uk/chembl/compound/inspect/",link,"', target = '_blank'>",database,"</a>")) %>% 
    select(internal_id, link, external_id, database) %>% 
    distinct()
  
  dgidb.links <- read.table(syn$get("syn12684108")$path, sep = "\t", quote = "", header = T) %>% 
    select(drug_claim_primary_name, drug_claim_name, drug_name, drug_chembl_id) %>% 
    distinct() %>% 
    mutate(drug_claim_primary_name2 = drug_claim_primary_name) %>% 
    gather("namesource", "name", -drug_claim_primary_name) %>% 
    select(-namesource) %>% 
    filter(!grepl("^\\d+$", name) & name != "") %>% 
    set_names(c("external_id", "common_name")) %>% 
    distinct() %>% 
    left_join(db.names) %>% 
    mutate(link = sapply(external_id, function(x) URLencode(x))) %>%
    mutate(link = paste0("<a href='http://www.dgidb.org/interaction_search_results/",link,"', target = '_blank'>",database,"</a>")) %>% 
    filter(!is.na(internal_id)) %>% 
    select(internal_id, link, external_id, database) %>% 
    distinct()
  
  drugbank.links <- db.names %>% 
    filter(database == "drugbank") %>% 
    mutate(link = sapply(external_id, function(x) URLencode(x))) %>%
    mutate(link = paste0("<a href='https://www.drugbank.ca/drugs/",link,"', target = '_blank'>",database,"</a>")) %>% 
    select(internal_id, link, external_id, database) %>% 
    distinct()
  
  db.links <- bind_rows(chembl.links, drugbank.links, dgidb.links)
  db.links$link <- gsub(">chembl<", ">ChEMBL<", db.links$link)
  db.links$link <- gsub(">dgidb<", ">DGIdb<", db.links$link)
  db.links$link <- gsub(">drugbank<", ">DrugBank<", db.links$link)
}