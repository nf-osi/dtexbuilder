synapseclient <- NULL

.onLoad <- function(libname, pkgname) {
  
  syn_inst <- reticulate::py_module_available("synapseclient")
  
  if(syn_inst){
    synapseclient <<- reticulate::import('synapseclient', delay_load = T)
  }
  
}