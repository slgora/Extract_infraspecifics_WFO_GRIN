# Tools/clean_outlink_tokens.R
# Filter semicolon-separated OUTLINK_PAGE tokens to keep only bona fide taxon-like strings.
# Functions:
#  - clean_outlink_field(outcell, species_binomial)
#  - clean_outlinks_in_df(df)

library(stringr)
library(dplyr)

is_taxon_token <- function(token, genus = NULL) {
  if (is.na(token) || !nzchar(token)) return(FALSE)
  tok <- str_squish(token)
  ui_blacklist <- c("\\band\\b", "\\bimplication", "\\ball names\\b", "\\baccepted names\\b",
                    "\\bclick here\\b", "\\bexportable\\b", "\\bexplore the data\\b", "\\bcheck a plant name\\b",
                    "\\bapplication of the\\b", "\\band its implications\\b", "\\band allied genera\\b")
  if (any(str_detect(tolower(tok), ui_blacklist))) return(FALSE)
  if (!is.null(genus) && nzchar(genus)) {
    genus_cap <- paste0(toupper(substr(genus,1,1)), tolower(substr(genus,2,nchar(genus))))
    if (!str_detect(tok, paste0("^(", genus_cap, "|", tolower(genus), ")\\b"))) return(FALSE)
  }
  pattern <- paste0(
    "^",
    "([A-Z][a-z]+|[a-z]+)\\s+[a-z\\-]+",
    "(?:\\s+(?:subsp\\.?|subspecies|subvar\\.?|var\\.?|variety|forma|f\\.?|form\\b|cv\\.?|cultivar|race)\\s+[a-z\\-]+|\\s+[a-z\\-]+)?",
    "(?:\\s+[A-Za-z\\.\\-\\&'()]+)*",
    "$"
  )
  ok <- str_detect(tok, regex(pattern))
  ok && length(str_split(tok, "\\s+")[[1]]) >= 2
}

clean_outlink_field <- function(outcell, species_binomial) {
  if (is.na(outcell) || !nzchar(outcell)) return(character(0))
  parts <- unlist(str_split(as.character(outcell), "\\s*;\\s*"))
  parts <- unique(trimws(parts))
  genus <- str_split(str_squish(species_binomial), "\\s+")[[1]][1]
  keep <- vapply(parts, function(p) is_taxon_token(p, genus = genus), logical(1))
  parts[keep]
}

clean_outlinks_in_df <- function(df) {
  if (!("input_species" %in% names(df)) || !("infraspecific_list_OUTLINK_PAGE" %in% names(df))) {
    stop("Data frame must contain 'input_species' and 'infraspecific_list_OUTLINK_PAGE' columns")
  }
  
  df2 <- df
  df2$OUTLINK_extracted_clean <- vector("list", nrow(df2))
  df2$n_infraspecific_OUTLINK_PAGE_clean <- integer(nrow(df2))
  
  for (i in seq_len(nrow(df2))) {
    sp <- as.character(df2$input_species[i])
    cell <- df2$infraspecific_list_OUTLINK_PAGE[i]
    cleaned_tokens <- clean_outlink_field(cell, sp)
    df2$OUTLINK_extracted_clean[[i]] <- cleaned_tokens
    df2$n_infraspecific_OUTLINK_PAGE_clean[i] <- length(cleaned_tokens)
  }
  
  canonicalize_infra <- function(name) {
    if (is.na(name) || !nzchar(name)) return(NA_character_)
    nm <- str_squish(name)
    toks <- str_split(nm, "\\s+")[[1]]
    toks_clean <- tolower(gsub("[[:punct:]]", "", toks))
    if (length(toks_clean) < 3) return(NA_character_)
    rank_tokens <- c("subsp","subspecies","subvar","var","variety","forma","form","f","cultivar","cv","race")
    rank_idx <- which(toks_clean %in% rank_tokens)
    if (length(rank_idx) >= 1) {
      idx <- rank_idx[1]
      if (length(toks_clean) >= idx + 1) return(paste(toks_clean[1], toks_clean[2], toks_clean[idx], toks_clean[idx + 1], sep = " "))
    }
    third <- toks[3]
    if (grepl("^[a-z]", third)) return(paste(toks_clean[1], toks_clean[2], toks_clean[3], sep = " "))
    NA_character_
  }
  
  df2$infraspecific_list_OUTLINK_PAGE_clean <- sapply(df2$OUTLINK_extracted_clean, function(x) if (length(x) == 0) NA_character_ else paste(sort(unique(x)), collapse = "; "))
  df2$n_infraspecific_OUTLINK_PAGE <- df2$n_infraspecific_OUTLINK_PAGE_clean
  
  df2$infraspecific_list_combined <- NA_character_
  df2$n_infraspecific_combined <- 0
  for (i in seq_len(nrow(df2))) {
    grin <- if ("infraspecific_list_GRIN" %in% names(df2)) df2$infraspecific_list_GRIN[i] else NA_character_
    wfo  <- if ("infraspecific_list_WFO" %in% names(df2)) df2$infraspecific_list_WFO[i] else NA_character_
    wfo_inc <- if ("infraspecific_list_WFO_included_variety" %in% names(df2)) df2$infraspecific_list_WFO_included_variety[i] else NA_character_
    wfo_syn <- if ("infraspecific_list_WFO_synonyms" %in% names(df2)) df2$infraspecific_list_WFO_synonyms[i] else NA_character_
    out_clean <- df2$OUTLINK_extracted_clean[[i]]
    
    split_and_keep <- function(cell) {
      if (is.na(cell) || !nzchar(cell)) return(character(0))
      parts <- unlist(str_split(as.character(cell), "\\s*;\\s*"))
      unique(trimws(parts))
    }
    
    grin_found <- split_and_keep(grin)
    wfo_found  <- split_and_keep(wfo)
    wfo_inc_found <- split_and_keep(wfo_inc)
    wfo_syn_found <- split_and_keep(wfo_syn)
    
    all_candidates <- unique(c(grin_found, wfo_found, wfo_inc_found, wfo_syn_found, out_clean))
    keys <- vapply(all_candidates, canonicalize_infra, FUN.VALUE = character(1))
    keep_idx <- which(!is.na(keys))
    if (length(keep_idx) == 0) {
      df2$infraspecific_list_combined[i] <- NA_character_
      df2$n_infraspecific_combined[i] <- 0
      next
    }
    keys_clean <- keys[keep_idx]
    candidates_clean <- all_candidates[keep_idx]
    unique_keys <- unique(keys_clean)
    key_to_preferred <- function(k) {
      idx <- which(vapply(wfo_inc_found, function(x) identical(canonicalize_infra(x), k), logical(1)))
      if (length(idx) > 0) return(wfo_inc_found[idx[1]])
      idx <- which(vapply(wfo_found, function(x) identical(canonicalize_infra(x), k), logical(1)))
      if (length(idx) > 0) return(wfo_found[idx[1]])
      idx <- which(vapply(grin_found, function(x) identical(canonicalize_infra(x), k), logical(1)))
      if (length(idx) > 0) return(grin_found[idx[1]])
      idx <- which(vapply(out_clean, function(x) identical(canonicalize_infra(x), k), logical(1)))
      if (length(idx) > 0) return(out_clean[idx[1]])
      idx <- which(vapply(candidates_clean, function(x) identical(canonicalize_infra(x), k), logical(1)))
      if (length(idx) > 0) return(candidates_clean[idx[1]])
      NA_character_
    }
    combined_pref <- unique(na.omit(vapply(unique_keys, key_to_preferred, FUN.VALUE = character(1), USE.NAMES = FALSE)))
    df2$infraspecific_list_combined[i] <- if (length(combined_pref) == 0) NA_character_ else paste(sort(combined_pref), collapse = "; ")
    df2$n_infraspecific_combined[i] <- length(combined_pref)
  }
  
  df2$OUTLINK_extracted_clean <- NULL
  df2$n_infraspecific_OUTLINK_PAGE_clean <- NULL
  df2
}
