#' Get an InChI value from an InChIKey from ChemSpider.
#' @description Uses the ChemSpider REST API, this service does not require an API key. I think it will always provide a standard InChI?
#' @param inchikey An InChIKey. 
#' @return character, InChI. 
#' @export
#'
.convert_inchikey_to_inchi_chemspider <- function(inchikey){
  Sys.sleep(0.25) ##to prevent requests from happening too fast
  
  statement <- glue::glue('www.chemspider.com/InChI.asmx/InChIKeyToInChI?inchi_key={inchikey}')
  
  res <- httr::GET(statement)
  
  if(res$status_code==200){
    inchi <- XML::xmlToList(rawToChar(res$content))
    if(is.null(inchi)){inchi <- NA}
  }else{
    message(glue::glue('inchikey "{inchikey}" appears to be invalid'))
    inchi <- NA
  }
  inchi
  
}

#' Get an InChI value from an InChIKey via Pubchem
#' @description Uses the PubChem PUG REST API, this service does not require an API key. I think it will always provide a standard InChI?
#' @param inchikey An InChIKey. 
#' @return character, InChI. 
#' @export
#'
.convert_inchikey_to_inchi_pubchem <- function(inchikey){
  Sys.sleep(0.25) ##to prevent requests from happening too fast
  
  statement <- glue::glue('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/InChI/xml')
  
  res <- httr::with_config(httr::config(http_version = 0), {
    httr::GET(statement)
    })
  
  if(res$status_code==200){
    inchi <- XML::xmlToList(rawToChar(res$content))$Properties$InChI
    if(is.null(inchi)){inchi <- NA}
  }else{
    message(glue::glue('inchikey "{inchikey}" appears to be invalid'))
    inchi <- NA
  }
  inchi

}

#' Get an InChI value from an InChIKey via Pubchem
#' @description Uses the PubChem PUG REST API, this service does not require an API key. I think it will always provide a standard InChI?
#' @param inchikey An InChIKey. 
#' @return character, InChI. 
#' @export
#'
.convert_id_to_structure_pubchem <- function(input_id, id_type = c("name", "inchikey"), output_type = c("InChI", "InChIKey", "CanonicalSMILES", "IsomericSMILES")){
  Sys.sleep(0.25) ##to prevent requests from happening too fast
  
  input <- URLencode(input_id)
  
  statement <- glue::glue('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{id_type}/{input}/property/{output_type}/xml')
  
  res <- httr::with_config(httr::config(http_version = 0), {
    httr::GET(statement)
  })
  
  if(res$status_code==200){
    res_2 <- XML::xmlToList(rawToChar(res$content))
    struct <- purrr::pluck(res_2, "Properties", 2)
    if(is.null(struct)){struct <- NA}
  }else{
    message(glue::glue('input "{input_id}" appears to be invalid'))
    struct <- NA
  }
  
  struct
  
}
