#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
process_structures <- function(...){
  
  dfs <- purrr::map(list(...), dplyr::mutate_all, as.character)
  
  structures <- dplyr::bind_rows(dfs)
  
  inchikey_smiles <- structures  %>% 
    dplyr::select(inchikey, inchi, smiles) %>% 
    dplyr::group_by(inchikey) %>% 
    dplyr::slice(1) %>% 
    dplyr::ungroup()
  
  #TIf there are more than one InChI or SMILES per InChIKey group (most likely from some of the hand-curated data sources), take the first InChI/SMILES in each group; this will result in a preference for the the 
  #ChEMBL-assigned InChI/SMILES, if it exists, which is what I prefer.
  
  message('validating smiles...')
  
  parser <- rcdk::get.smiles.parser()
  valid.smiles<- pbapply::pbsapply(inchikey_smiles$smiles, webchem::is.smiles)
  
  valid.smiles <- as.data.frame(valid.smiles)
  valid.smiles$smiles <- inchikey_smiles$smiles
  valid.smiles$inchikey <- inchikey_smiles$inchikey
  valid.smiles <- dplyr::distinct(valid.smiles)
  
  message('standardizing valid smiles...')
  
  valid <- valid.smiles$smiles[valid.smiles$valid.smiles==TRUE]
  
  if(!reticulate::py_module_available('molvs')){
    message("molvs python module not available. install using reticulate::py_install('molvs')")
  }else{
  
  molvs <- reticulate::import('molvs')
  
  standardized_smiles <- pbapply::pbsapply(unique(valid), function(x){
    tryCatch({
      molvs$standardize_smiles(reticulate::r_to_py(x))
    }, warning = function(w) {
      message(paste("structure",x,"has an issue"))
    }, error = function(e) {
      message(paste("structure",x,"is invalid"))
    })
  })
  
  names(standardized_smiles) <- unique(valid)
  
  valid.df <- standardized_smiles %>% 
    plyr::ldply() %>%
    purrr::set_names(c("smiles", "std_smiles")) %>% 
    dplyr::filter(!is.na(std_smiles)) %>% 
    dplyr::inner_join(valid.smiles) %>% 
    dplyr::select(smiles, std_smiles, inchikey)
  
  structure_ids <- dplyr::select(structures, database, external_id, inchi, inchikey) %>% 
    dplyr::distinct()
    
  structures_with_std_smiles <- dplyr::inner_join(structure_ids, valid.df)
  
  }  
}


#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
.parse_fingerprint <- function(input, type) {
  print("parsing smiles")
  input.mol <- rcdk::parse.smiles(as.character(input))
  print("doing typing")
  pbapply::pblapply(input.mol, rcdk::do.typing)
  print("doing aromaticity")
  pbapply::pblapply(input.mol, rcdk::do.aromaticity)
  print("doing isotopes")
  pbapply::pblapply(input.mol, rcdk::do.isotopes)
  print("generating fingerprints")
  pbapply::pblapply(input.mol, rcdk::get.fingerprint, type = type)
}

#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
generate_fingerprints <- function(structure_df, type){ 
  
  valid <- as.character(structure_df$std_smiles)
  
  parser <- rcdk::get.smiles.parser()
  foo <- list()
  ct <- 1

for(i in 1:ceiling(length(valid)/5000)){
  if((length(valid)-(i*5000))>=0){
    print(ct)
    print(i*5000)
    print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
    foo <- append(foo, .parse_fingerprint(valid[ct:(i*5000)], type = type))
    ct<-ct+5000
  }else{
    print(ct)
    print(length(valid))
    print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
    foo <- append(foo, .parse_fingerprint(valid[ct:length(valid)], type = type))
  }
}
  names(foo) <- structure_df$inchikey
  foo
}
  
#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
process_activities <- function(processed_structures, ...){
  
  ext_id_inchikey <- dplyr::select(processed_structures, external_id, inchikey) %>% 
    dplyr::distinct()
  
  dfs <- purrr::map(list(...), function(x){
    x %>% 
      dplyr::mutate_all(as.character)
  })
  
  activities <- dplyr::bind_rows(dfs) %>% 
    dplyr::mutate(dplyr::across(c('pchembl_value','standard_value'), as.numeric))
  
  
  message("processing quantitative values...")
  ##quantitative data
  quant <- activities %>% dplyr::filter(evidence_type == "quantitative")
  message("summarizing pchembl values...")
  
  pchembl_summary <- quant %>% 
    dplyr::left_join(x=ext_id_inchikey, y = .) %>% 
    dplyr::select(inchikey, hugo_gene, pchembl_value) %>% 
    dplyr::filter(!is.na(pchembl_value)) %>% 
    dplyr::filter(pchembl_value != "NULL") %>% 
    dplyr::group_by(inchikey, hugo_gene) %>% 
    dplyr::summarize("n_quantitative" = dplyr::n(), 
              "mean_pchembl" = mean(pchembl_value, na.rm = T),
              "cv" = raster::cv(pchembl_value),
              "sd" = sd(pchembl_value)) %>% 
    dplyr::ungroup()
  
  message("summarizing assay values...")
  
  chembl_assaytype_summary <- quant %>%
    dplyr::left_join(x=ext_id_inchikey, y = .) %>% 
    dplyr::select(inchikey, hugo_gene, standard_value, standard_type) %>% 
    dplyr::filter(standard_type %in% c("IC50", "AC50", "EC50", "C50", "Potency", "Ki", "Kd")) %>% 
    dplyr::filter(!is.na(standard_type)) %>% 
    dplyr::group_by(inchikey, hugo_gene, standard_type) %>% 
    dplyr::summarize(mean_value = mean(standard_value)) %>%
    dplyr::ungroup() %>% 
    tidyr::spread(standard_type, mean_value) %>% 
    dplyr::select(inchikey, hugo_gene, IC50, AC50, EC50, Potency, Ki, Kd) %>% 
    purrr::set_names(c("inchikey", "hugo_gene", "IC50_nM","AC50_nM",
                "EC50_nM", "Potency_nM", "Ki_nM", "Kd_nM"))
  
  quant_data <- dplyr::full_join(pchembl_summary, chembl_assaytype_summary)
  
  ##qualitative data
  
  message("summarizing qualitative associations...")
  
  qual <- activities %>% dplyr::filter(evidence_type == "qualitative")
  

  qual_data <- qual %>% 
    dplyr::left_join(x=ext_id_inchikey, y = .) %>% 
    dplyr::select(inchikey, hugo_gene) %>% 
    dplyr::group_by(inchikey, hugo_gene) %>% 
    dplyr::count(name = "n_qualitative") %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(!is.na(hugo_gene)) %>% 
    dplyr::distinct()
    
  full.db <- dplyr::full_join(quant_data, qual_data, by = c("hugo_gene", "inchikey"))
  
  ####need to add confidence calculations...
  
  message('calculating confidence and selectivity scores...')
  
  total_qualitative <- sum(full.db$n_qualitative, na.rm = T)
  total_quantitative <- sum(full.db$n_quantitative, na.rm = T)
  
  full.db.scored <- full.db %>% 
    dplyr::group_by(inchikey, hugo_gene) %>% 
    dplyr::mutate(total_n = sum(n_quantitative,n_qualitative, na.rm=T)) %>% 
    dplyr::ungroup()
  
  sd.n <- sd(full.db.scored$total_n, na.rm = T)
  mean.n <- mean(full.db.scored$total_n, na.rm = T)
  
  full.db.scored <- full.db.scored %>% 
    dplyr::mutate(confidence = (total_n-mean.n)/sd.n) %>% 
    dplyr::group_by(inchikey) %>% 
    dplyr::mutate(pchembl_d = sum(mean_pchembl, na.rm = T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(hugo_gene) %>% 
    dplyr::mutate(pchembl_t = sum(mean_pchembl, na.rm = T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(inchikey, hugo_gene) %>% 
    # mutate(ksi_dt = mean_pchembl/((pchembl_t+pchembl_d)-mean_pchembl)) %>%  ##normalizes for target frequency, perhaps not useful for "drug" selectivity
    dplyr::mutate(known_selectivity_index = mean_pchembl/pchembl_d) %>% 
    dplyr::mutate(confidence = signif(confidence,3)) %>% 
    dplyr::mutate(mean_pchembl = signif(mean_pchembl,3)) %>% 
    dplyr::mutate(known_selectivity_index = signif(known_selectivity_index,3)) %>% 
    dplyr::ungroup()
}

#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
process_names <- function(processed_structures, processed_activities, ...){
  
  ext_id_inchikey <- dplyr::select(processed_structures, external_id, inchikey) %>% 
    dplyr::distinct() %>% 
    dplyr::filter(inchikey %in% processed_activities$inchikey)
  
  dfs <- purrr::map(list(...), function(x){
    x %>% 
      dplyr::mutate_all(as.character)
  })
  
  names <- dplyr::bind_rows(dfs) %>% 
    dplyr::inner_join(ext_id_inchikey)
  
  synonyms <- dplyr::select(names, inchikey, cmpd_name) %>% 
    dplyr::distinct() %>%
    dplyr::group_by(inchikey) %>% 
    dplyr::mutate(pref_name = cmpd_name[1]) %>% ##first occurrence of name is the "preferred name" 
    dplyr::mutate(upper=toupper(cmpd_name),
                  dupl=duplicated(upper,fromLast=FALSE)) %>%
    dplyr::filter(dupl!=TRUE) %>%
    dplyr::select(-upper,-dupl) %>% 
    dplyr::rename(synonym = cmpd_name) %>% 
    dplyr::ungroup()
  
  distinct_names <- synonyms %>% 
    dplyr::select(inchikey, pref_name) %>% 
    dplyr::distinct()
  
  list("names" = distinct_names, "synonyms" = synonyms)
}

#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
filter_structures <- function(processed_structures, processed_activities){
  
  processed_structures %>% 
    dplyr::select(inchikey, inchi, std_smiles) %>% 
    dplyr::filter(inchikey %in% processed_activities$inchikey)  %>% 
    dplyr::distinct()
  
}

#' Placeholder
#' @description Placeholder
#' @param placeholder
#' @return TBD
#' @export
#' 
process_external_ids <- function(processed_structures, processed_activities){
  
  ext_id_inchikey <- dplyr::select(processed_structures, external_id, inchikey) %>% 
    dplyr::filter(inchikey %in% processed_activities$inchikey) %>% 
    dplyr::distinct()
  
}


