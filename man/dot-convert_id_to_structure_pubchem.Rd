% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conversion_utils.R
\name{.convert_id_to_structure_pubchem}
\alias{.convert_id_to_structure_pubchem}
\title{Get an InChI value from an InChIKey via Pubchem}
\usage{
.convert_id_to_structure_pubchem(
  input_id,
  id_type = c("name", "inchikey"),
  output_type = c("InChI", "InChIKey", "CanonicalSMILES", "IsomericSMILES")
)
}
\arguments{
\item{inchikey}{An InChIKey.}
}
\value{
character, InChI.
}
\description{
Uses the PubChem PUG REST API, this service does not require an API key. I think it will always provide a standard InChI?
}
