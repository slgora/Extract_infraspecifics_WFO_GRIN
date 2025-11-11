# Functions/Query_taxa_resolver.R
# Query Global Names Verifier API for a single taxon.
# Uses all_matches=true & species_group=true to surface infraspecific candidates.
# Returns parsed JSON (list) or NULL on error.

library(httr)
library(jsonlite)

query_taxa_resolver <- function(taxa, sources = c('196'), verbose = FALSE, timeout_sec = 30) {
  if (!is.character(taxa) || length(taxa) != 1) {
    if (verbose) message("query_taxa_resolver: 'taxa' must be a single character string")
    return(NULL)
  }
  
  taxa_encoded <- URLencode(taxa, reserved = TRUE)
  URL <- paste0(
    'https://verifier.globalnames.org/api/v1/verifications/', taxa_encoded,
    '?data_sources=', paste(sources, collapse = "|"),
    '&all_matches=true&capitalize=true&species_group=true',
    '&fuzzy_uninomial=false&stats=false&main_taxon_threshold=0.8'
  )
  
  if (verbose) message("GET ", URL)
  res <- try(GET(URL, user_agent("R (taxa standardization script)"), timeout(timeout_sec)), silent = TRUE)
  if (inherits(res, "try-error")) {
    if (verbose) message("query_taxa_resolver: GET error - ", conditionMessage(attr(res, "condition")))
    return(NULL)
  }
  if (status_code(res) != 200) {
    if (verbose) message("query_taxa_resolver: HTTP status ", status_code(res))
    return(NULL)
  }
  txt <- content(res, "text", encoding = "UTF-8")
  parsed <- try(fromJSON(txt, simplifyVector = FALSE), silent = TRUE)
  if (inherits(parsed, "try-error")) {
    if (verbose) message("query_taxa_resolver: JSON parse error")
    return(NULL)
  }
  parsed
}
