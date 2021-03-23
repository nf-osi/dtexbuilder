#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
process_drugbank <- function(db_synapse_id, structures_synapse_id){
  
  drugbank <- .get_drugbank_data(db_synapse_id, structures_synapse_id)
  
  message('processing activity data...')
  act <- .format_drugbank_activities(drugbank$db)
  
  message('processing structure data...')
  struct <- .format_drugbank_structures(drugbank$structures)
  
  message('processing names data...')
  names <- .format_drugbank_names(drugbank$structures)
  
  list("activities" = act, "structures" = struct, "names" = names)
  
}

#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
.format_drugbank_structures <- function(drugbank_structures){
  
  df <- drugbank_structures %>%
    dplyr::mutate(external_id = glue::glue('drugbank:{`DrugBank ID`}')) %>% 
    dplyr::rename(cmpd_name = Name,
                  inchikey = InChIKey,
                  inchi = InChI,
                  smiles = SMILES) %>% 
    dplyr::select(external_id, cmpd_name, inchikey, inchi, smiles) %>% 
    dplyr::mutate(database = "drugbank") %>% 
    dplyr::distinct() %>% 
    dplyr::filter(!is.na(inchikey) & !is.na(smiles))
    
  df 
  
}

#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
.format_drugbank_activities <- function(drugbank_db, drugbank_structures){
  
  message("retrieving activities and targets from database...")
  drugbank.summary <- lapply(1:(length(drugbank_db)-1), function(x){
    targets <- c()
    for(i in 1:length(drugbank_db[[x]]$targets)){
      targets<-c(targets, drugbank_db[[x]]$targets[[i]]$polypeptide$`external-identifiers`$`external-identifier`$identifier)
    }
    
    external_id <- c(rep(drugbank_db[[x]]$`drugbank-id`$text, length(targets)))
    
    data.frame(external_id, targets)
  })
  
  drugbank.df <- plyr::ldply(drugbank.summary) 
  
  message('querying ensembl for HGNC ids...')
  
  ensembl <- biomaRt::useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
  
  map <- biomaRt::getBM(attributes = c("hgnc_id", "hgnc_symbol"), filters = "hgnc_id", values = unique(drugbank.df$targets), mart = ensembl) %>% 
    purrr::set_names(c("targets", "hugo_gene"))
  
  drugbank.df <- dplyr::inner_join(drugbank.df, map) %>% 
    dplyr::select(external_id, hugo_gene) %>% 
    tibble::add_column(evidence_type = "qualitative") 
  
}

#' Export the DrugBank names.
#' @description Export the DrugBank names.
#' @param placeholder
#' @return TBD
#' @export
#' 
.format_drugbank_names <- function(drugbank_structures){
  
  df <- drugbank_structures %>%
    dplyr::mutate(external_id = glue::glue('drugbank:{`DrugBank ID`}')) %>% 
    dplyr::rename(cmpd_name = Name) %>% 
    dplyr::select(external_id, cmpd_name) %>% 
    dplyr::distinct()
  
  df 
  
}

#' Retrieve drugbank data from database. 
#' @description Retrieve drugbank data from exported XML and csvs. Files should be hosted on synapse. 
#' @param db_synapse_id The Synapse ID of the drugbank XML. 
#' @param structures_synapse_id The Synapse ID of the structures CSV. 
#' @return TBD
#' @export
#' 
.get_drugbank_data <- function(db_synapse_id, structures_synapse_id){
  
  .check_login()
  
  message('reading xml database - expect this to take a while...')
  db <- .syn$get(db_synapse_id)$path %>% 
    XML::xmlToList()
  
  message('reading structures file...')
  structures <- .syn$get(structures_synapse_id)$path %>% 
    readr::read_csv()
  
  list("db" = db, "structures" = structures)
}





