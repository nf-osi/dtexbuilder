#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
generate_external_links <- function(external_ids, activities){
  
  cl <- .cmpd_links(external_ids = external_ids)
  
  gl <- .gene_links(activities = activities)
  
  list('cmpd_links' = cl, 'gene_links' = gl)
}



.cmpd_links <- function(external_ids){

  chembl_link <- "<a href='https://www.ebi.ac.uk/chembl/compound/inspect/{link_string}', target = '_blank'>ChEMBL</a>"
  dgidb_link <- "<a href='http://www.dgidb.org/interaction_search_results/{link_string}', target = '_blank'>DGIdb</a>"
  drugbank_link <- "<a href='https://go.drugbank.com/drugs/{link_string}', target = '_blank'>DrugBank</a>"
  chemicalprobes_link <- "<a href='https://www.chemicalprobes.org/{link_string}', target = '_blank'>ChemicalProbes</a>"
    
  links <- external_ids %>% 
    dplyr::select(external_id, inchikey, database) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(link_string = stringr::str_remove(external_id, ".+?\\:")) %>% 
    dplyr::mutate(link_string = sapply(link_string, function(x) URLencode(x))) %>%
    dplyr::mutate(link = dplyr::case_when(database == "chembl" ~ glue::glue(chembl_link),
                                   database == "dgidb" ~ glue::glue(dgidb_link),
                                   database == "chemicalprobes" ~ glue::glue(chemicalprobes_link),
                                   database == "drugbank" ~ glue::glue(drugbank_link))) %>% 
    dplyr::select(-link_string)
  
}


.gene_links <- function(activities){
  
  gene_link <- "<a href = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene={link_string}', target = '_blank'>GeneCards</a>"
    
  genes <- unique(activities$hugo_gene) %>% 
    dplyr::as_tibble() %>% 
    purrr::set_names(c("hugo_gene")) %>% 
    dplyr::mutate(link_string = sapply(hugo_gene, function(x) URLencode(x))) %>%
    dplyr::mutate(link = glue::glue(gene_link)) %>% 
    dplyr::select(-link_string)
  
}