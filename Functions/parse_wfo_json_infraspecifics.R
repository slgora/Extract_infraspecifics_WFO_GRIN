# Functions/parse_wfo_json_infraspecifics.R
# Parse WFO RDF/JSON payloads returned by list.worldfloraonline.org outlinks.
# Extracts:
#  - included_varieties : nodes referenced by hasPart/haspart-like predicates
#  - synonyms           : nodes referenced by hasSynonym/isSynonym-like predicates
#  - page_candidates    : any fullName literals that look infraspecific

library(jsonlite)
library(stringr)

.get_prop_values <- function(node, pattern) {
  if (is.null(node) || !is.list(node)) return(character(0))
  key_names <- names(node); if (is.null(key_names)) return(character(0))
  hits <- key_names[str_detect(tolower(key_names), tolower(pattern))]
  vals <- character(0)
  for (k in hits) {
    v <- node[[k]]
    if (is.list(v)) {
      for (el in v) if (is.list(el) && !is.null(el$value)) vals <- c(vals, as.character(el$value))
    } else if (is.character(v)) vals <- c(vals, v)
  }
  unique(vals)
}

.resolve_uri_to_node <- function(parsed_json, uri) {
  if (is.null(parsed_json) || !nzchar(uri)) return(NULL)
  if (!is.null(parsed_json[[uri]])) return(parsed_json[[uri]])
  keys <- names(parsed_json); if (is.null(keys)) return(NULL)
  simple <- sub("^.*/", "", uri)
  k <- keys[str_detect(keys, fixed(simple))]
  if (length(k) > 0) return(parsed_json[[k[1]]])
  NULL
}

.is_infraspecific_name <- function(name, species_binomial = NULL) {
  if (is.null(name) || !nzchar(name)) return(FALSE)
  nm <- str_squish(name); toks <- str_split(nm, "\\s+")[[1]]
  if (length(toks) < 3) return(FALSE)
  if (!is.null(species_binomial) && nzchar(species_binomial)) {
    bin <- str_split(str_squish(species_binomial), "\\s+")[[1]]
    if (length(bin) >= 2 && !(tolower(toks[1]) == tolower(bin[1]) && tolower(toks[2]) == tolower(bin[2]))) return(FALSE)
  }
  third <- toks[3]
  if (grepl("^(?i)(subsp\\.?|subspecies|subvar\\.?|var\\.?|variety|forma|f\\.?|form\\b|cv\\.?|cultivar|race)$", third, perl = TRUE)) return(TRUE)
  if (grepl("^[a-z]", third)) return(TRUE)
  FALSE
}

parse_wfo_json_infraspecifics <- function(raw_json_text, species_binomial = NULL, verbose = FALSE) {
  if (is.null(raw_json_text) || !nzchar(raw_json_text)) return(list(included_varieties = character(0), synonyms = character(0), page_candidates = character(0)))
  parsed <- try(fromJSON(raw_json_text, simplifyVector = FALSE), silent = TRUE)
  if (inherits(parsed, "try-error")) return(list(included_varieties = character(0), synonyms = character(0), page_candidates = character(0)))
  keys <- names(parsed); if (is.null(keys) || length(keys) == 0) return(list(included_varieties = character(0), synonyms = character(0), page_candidates = character(0)))
  
  included <- character(0); syns <- character(0); candidates <- character(0)
  
  for (k in keys) {
    node <- parsed[[k]]
    fn <- .get_prop_values(node, "fullname|full name|fullName")
    if (length(fn) > 0) for (f in fn) if (.is_infraspecific_name(f, species_binomial)) candidates <- c(candidates, f)
  }
  
  for (k in keys) {
    node <- parsed[[k]]
    parts <- .get_prop_values(node, "haspart|hasPart|hasMember")
    if (length(parts) > 0) {
      for (p in parts) {
        ref <- .resolve_uri_to_node(parsed, p)
        if (!is.null(ref)) {
          fn2 <- .get_prop_values(ref, "fullname|full name|fullName")
          if (length(fn2) > 0) for (f in fn2) if (.is_infraspecific_name(f, species_binomial)) included <- c(included, f)
        }
      }
    }
  }
  
  for (k in keys) {
    node <- parsed[[k]]
    syn_uris <- .get_prop_values(node, "hassynonym|hasSynonym|isSynonym|has_synonym")
    if (length(syn_uris) > 0) {
      for (su in syn_uris) {
        ref <- .resolve_uri_to_node(parsed, su)
        if (!is.null(ref)) {
          fn2 <- .get_prop_values(ref, "fullname|full name|fullName")
          if (length(fn2) > 0) for (f in fn2) if (.is_infraspecific_name(f, species_binomial)) syns <- c(syns, f)
        }
      }
    }
  }
  
  list(included_varieties = unique(included), synonyms = unique(syns), page_candidates = unique(candidates))
}
