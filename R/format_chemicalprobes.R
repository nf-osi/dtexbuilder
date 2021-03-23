#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
process_chemicalprobes <- function(synapse_id){
  
  chemicalprobes_df <- .get_chemicalprobes_data(synapse_id)
  
  message('processing activity data...')
  act <- .format_chemicalprobes_activities(chemicalprobes_df)

  message('processing structure data...')
  struct <- .format_chemicalprobes_structures(chemicalprobes_df)
  
  message('processing names data...')
  names <- .format_chemicalprobes_names(chemicalprobes_df)
  
  list("activities" = act, "structures" = struct, "names" = names)
  
}


#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
.format_chemicalprobes_structures <- function(chemicalprobes_df){
  
  df <- chemicalprobes_df %>%
    dplyr::select(external_id, cmpd_name, inchikey, smiles) %>% 
    dplyr::mutate(database = "chemicalprobes") %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(inchi = sapply(inchikey, .convert_id_to_structure_pubchem, id_type = "inchikey", output_type = "InChI"))
  
  df_2 <- df %>% 
    dplyr::filter(is.na(inchi)) %>% 
    dplyr::mutate(inchi = sapply(inchikey, .convert_inchikey_to_inchi_chemspider))
  
  if(nrow(df_2) > 0){
    df_3 <- df %>% 
      dplyr::bind_rows(df_2) %>%
      dplyr::filter(!is.na(inchi)) %>% 
      dplyr::filter(!is.na(inchikey))
  }else{
    df_3 <- df
  }
  
  df_3 
  
}

#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
.format_chemicalprobes_activities <- function(chemicalprobes_df){
  
  chemicalprobes_df %>% 
    dplyr::select(external_id, cmpd_name, hugo_gene) %>% 
    tidyr::separate_rows(hugo_gene, sep = ",") %>% 
    dplyr::mutate(hugo_gene = stringr::str_trim(hugo_gene, side = "both")) %>% #catch and remove whitespace
    dplyr::distinct() %>% 
    tibble::add_column(evidence_type = "qualitative") 
  
}

#' Export the chemicalprobes names.
#' @description Export the chemicalprobes names.
#' @param placeholder
#' @return TBD
#' @export
#' 
.format_chemicalprobes_names <- function(chemicalprobes_df){
  
  chemicalprobes_df %>% 
    dplyr::select(external_id, cmpd_name) %>% 
    dplyr::distinct()
  
}

#' Retrieve chemicalprobes data from CSVs. 
#' @description Retrieve chemicalprobes data from CSV (db exports must be requested by email as of early 2021). Files should be hosted on synapse. 
#' @param synapse_id The Synapse ID of the ChemicalProbes.org exported csv. 
#' @return TBD
#' @export
#' 
.get_chemicalprobes_data <- function(synapse_id){
  
  .check_login()
  
  .syn$get(synapse_id)$path %>%
    readr::read_csv(comment = "") %>% 
    dplyr::mutate(external_id = glue::glue('chemicalprobes:{`Chemical Probes Portal ID`}'), .keep = 'unused') %>% 
    dplyr::rename(cmpd_name = `Probe Name`, inchikey = `InChi Key`, smiles = `SMILES string`, hugo_gene = `Target name`) 

}




