#' Organize chembl dataframes.
#' @description Organize chembl data and export in DTEX-mergeable format
#' @param activity_df
#' @param names_df
#' @param structures_df
#' @return TBD
#' @export
#' 
process_chembl <- function(activity_file_id, names_file_id, structures_file_id){
  
  message('getting files...')
  chembl <- .get_chembl_files(activity_file_id, names_file_id, structures_file_id)
  
  message('processing activity data...')
  act <- .format_chembl_activities(activity_df = chembl$activity, names_df = chembl$names)
  
  message('processing structure data...')
  struct <- .format_chembl_structures(activity_df = chembl$activity, structures_df = chembl$structures, names_df = chembl$names)
  
  message('processing names data...')
  names <- .format_chembl_names(names_df = chembl$names)
  
  list("activities" = act, "structures" = struct, "names" = names)
  
}


#' Organize chembl structures
#' @description Organize chembl data and export in DTEX-mergeable format
#' @param activity_df
#' @param structures_df
#' @return TBD
#' @export
#' 
.format_chembl_structures <- function(activity_df, structures_df, names_df){
  
  chembl_ids <- dplyr::select(names_df, molregno, chembl_id) %>% 
    dplyr::distinct()
  
  chembl_struct <- structures_df %>%
    dplyr::distinct() %>% 
    dplyr::filter(molregno %in% activity_df$molregno) %>% 
    dplyr::inner_join(chembl_ids) %>% 
    dplyr::mutate(external_id = glue::glue("chembl:{chembl_id}")) %>% 
    dplyr::rename(inchi = standard_inchi,
           inchikey = standard_inchi_key,
           smiles = canonical_smiles) %>% 
    dplyr::mutate(database = "chembl") 
  
}

#' Organize chembl dataframes.
#' @description Organize chembl data and export in DTEX-mergeable format
#' @param activity_df
#' @param names_df
#' @param structures_df
#' @return TBD
#' @export
#' 
.format_chembl_activities <- function(activity_df, names_df){
  
  chembl_ids <- dplyr::select(names_df, molregno, chembl_id) %>% 
    dplyr::distinct()
  
  ##activities
  chembl_activities <- activity_df %>% 
    dplyr::filter(organism == "Homo sapiens") %>% 
    dplyr::filter(pchembl_value != "NULL") %>% 
    dplyr::mutate(external_id = glue::glue("chembl:{chembl_id}")) %>% 
    dplyr::inner_join(chembl_ids) %>%
    dplyr::rename(hugo_gene = component_synonym) %>% 
    dplyr::mutate(evidence_type = "quantitative") %>% 
    dplyr::mutate(database = "chembl") %>% 
    dplyr::select(external_id, hugo_gene, pchembl_value, standard_relation, standard_value, standard_units, standard_type, evidence_type)
  
}

#' Organize chembl dataframes.
#' @description Organize chembl data and export in DTEX-mergeable format
#' @param activity_df
#' @param names_df
#' @param structures_df
#' @return TBD
#' @export
#' 
.format_chembl_names <- function(names_df){
  
    names <- names_df %>% 
    dplyr::mutate(external_id = glue::glue("chembl:{chembl_id}")) %>% 
    tidyr::pivot_longer(c(pref_name, synonyms, chembl_id), names_to = "colnames", values_to = "cmpd_name") %>% 
    dplyr::select(external_id, cmpd_name) %>% 
    dplyr::filter(cmpd_name != "NULL") %>% 
    dplyr::distinct()
  
}


#' Retrieve chembl data from TSVs. 
#' @description Retrieve chembl data from TSVs extracted from the mysql database. Files should be hosted on synapse. 
#' @param activity_file_id
#' @param names_file_id
#' @param structures_file_id
#' @return TBD
#' @export
#' 
.get_chembl_files <- function(activity_file_id, names_file_id, structures_file_id){
  
  .check_login()
  
  message("retrieving activity file...")
  activity <- .syn$get(activity_file_id)$path %>% 
    readr::read_tsv() %>% 
    dplyr::rename(assay_chembl_id = chembl_id)
  
  message("retrieving names file...")
  names <- .syn$get(names_file_id)$path %>% 
    readr::read_tsv()
  
  message("retrieving structures file...")
  structures <- .syn$get(structures_file_id)$path %>%
    readr::read_tsv(comment = "")

  dat <- list("activity" = activity, "names" = names, "structures" = structures)
  message("done!")
  
  dat
}
