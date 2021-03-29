#' Prepare dgidb drugs for structure mapping
#' @description Generates a formatted CSV that can be used as input into the PubChem Identifier Exchange Service (Input ID list = synonyms, operator type = same CID, Output IDs = InChIs,InChiKey,smiles)
#' @param dgidb_interactions_url A url of the drug gene interaction database, eg "https://www.dgidb.org/data/monthly_tsvs/2021-Jan/interactions.tsv"
#' @return An input tsv for the PubChem Identifier Exchange Service 
#' @export
#' 
prepare_dgidb_data <- function(dgidb_interactions_url, omit_sources = c("DrugBank", "ChemblInteractions")){
  
  readr::read_tsv(dgidb_interactions_url) %>% 
    dplyr::filter(!interaction_claim_source %in% omit_sources) %>% 
    dplyr::filter(!is.na(gene_name)) %>% 
    dplyr::select(drug_claim_primary_name, gene_name) %>% 
    purrr::set_names(c("cmpd_name_orig", "hugo_gene")) %>% 
    dplyr::mutate(database = "dgidb") 
}

#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
prepare_dgidb_names_for_pubchem <- function(dgidb_df){
  
  dgidb_df %>% 
    dplyr::select(cmpd_name_orig) %>% 
    dplyr::distinct() %>% 
    readr::write_tsv("dgidb_cmpd_names.tsv")
  
}



#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
process_dgidb <- function(dgidb_df, path_to_inchi_tsv, path_to_inchikey_tsv, path_to_smiles_tsv){
  
  message('processing structure data...')
  struct <- .format_dgidb_structures(dgidb_df = dgidb_df, path_to_inchi_tsv, path_to_inchikey_tsv, path_to_smiles_tsv)
  
  message('processing names data...')
  names <- .format_dgidb_names(dgidb_structures = struct)
  
  message('processing activity data...')
  act <- .format_dgidb_activities(dgidb_df = dgidb_df, dgidb_structures = struct)
  
  list("activities" = act, "structures" = struct, "names" = names)
  
}

#' This is a a draft function. 
#' @description Placeholder
#' @param placeholder
#' @return TBD
# .query_and_format_dgidb_structures <- function(dgidb_df){
#   
#   external_ids <- dgidb_df %>% 
#     dplyr::select(external_id, cmpd_name, database) %>% 
#     dplyr::distinct() 
#   
#   inchi <- external_ids %>% 
#     dplyr::mutate(inchi = sapply(cmpd_name, .convert_id_to_structure_pubchem, id_type = "name", output_type = "InChI"))
#   
#   inchikey <- readr::read_tsv(path_to_inchikey_tsv) %>% 
#     purrr::set_names(c("cmpd_name", "inchikey")) %>%
#     dplyr::filter(!is.na(inchikey)) %>% 
#     dplyr::group_by(cmpd_name) %>% 
#     dplyr::slice(1) %>% #some compounds have multiple, pick first occurrence as "real" InChIKey. This will certainly not be error free. 
#     dplyr::ungroup()
#   
#   smiles <- readr::read_tsv(path_to_smiles_tsv) %>% 
#     purrr::set_names(c("cmpd_name", "smiles")) %>%
#     dplyr::filter(!is.na(smiles)) %>% 
#     dplyr::group_by(cmpd_name) %>% 
#     dplyr::slice(1) %>% #some compounds have multiple, pick first occurrence as "real" smiles This will certainly not be error free. 
#     dplyr::ungroup()     
#   
#   df <- dplyr::inner_join(inchi, inchikey) %>% 
#     dplyr::inner_join(smiles) %>% 
#     dplyr::inner_join(external_ids)
# }


#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
.format_dgidb_structures <- function(dgidb_df, path_to_inchi_tsv, path_to_inchikey_tsv, path_to_smiles_tsv){

  inchi <- readr::read_tsv(path_to_inchi_tsv) %>% 
    purrr::set_names(c("cmpd_name_orig", "inchi")) %>%
    dplyr::filter(!is.na(inchi)) %>% 
    dplyr::group_by(cmpd_name_orig) %>% 
    dplyr::slice(1) %>% #some compounds have multiple, pick first occurrence as "real" InChI. This will certainly not be error free. 
    dplyr::ungroup()
  
  inchikey <- readr::read_tsv(path_to_inchikey_tsv) %>% 
    purrr::set_names(c("cmpd_name_orig", "inchikey")) %>%
    dplyr::filter(!is.na(inchikey)) %>% 
    dplyr::group_by(cmpd_name_orig) %>% 
    dplyr::slice(1) %>% #some compounds have multiple, pick first occurrence as "real" InChIKey. This will certainly not be error free. 
    dplyr::ungroup()
  
  smiles <- readr::read_tsv(path_to_smiles_tsv) %>% 
    purrr::set_names(c("cmpd_name_orig", "smiles")) %>%
    dplyr::filter(!is.na(smiles)) %>% 
    dplyr::group_by(cmpd_name_orig) %>% 
    dplyr::slice(1) %>% #some compounds have multiple, pick first occurrence as "real" smiles This will certainly not be error free. 
    dplyr::ungroup()     
  
  df <- dplyr::inner_join(inchi, inchikey) %>% 
    dplyr::inner_join(smiles) %>%
    dplyr::group_by(inchikey) %>% 
    dplyr::mutate(cmpd_name = cmpd_name_orig[1]) %>%  ##first occurrence of name is the "preferred name" 
    dplyr::mutate(external_id  = glue::glue("dgidb:{cmpd_name}")) %>% 
    dplyr::distinct() %>%
    dplyr::ungroup()
    
}

#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
.format_dgidb_activities <- function(dgidb_df, dgidb_structures){
  
  external_ids <- dgidb_structures %>% 
    dplyr::select(cmpd_name, external_id, cmpd_name_orig) %>% 
    dplyr::distinct()
  
  dgidb_df %>% 
    dplyr::inner_join(external_ids) %>% 
    tibble::add_column(evidence_type = "qualitative") %>% 
    dplyr::select(-cmpd_name_orig)
  
}

#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
.format_dgidb_names <- function(dgidb_structures){
  
  external_ids <- dgidb_structures %>% 
    dplyr::select(cmpd_name, external_id) %>% 
    dplyr::distinct()
  
}
