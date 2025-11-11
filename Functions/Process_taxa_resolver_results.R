# Functions/Process_taxa_resolver_results.R
# Extractor for infraspecific matches across verifier JSON shapes.
# Primary function:
#   extract_infraspecifics_single(res, input_species, source_name = NA,
#                                 accepted_only = TRUE, include_synonyms = FALSE, verbose = FALSE)
#
# Returns a data.frame/tibble with columns:
#   input_species, source, candidate_name, rank, status, matched_name, match_type, outlink

library(stringr)
library(dplyr)
library(tibble)

extract_infraspecifics_single <- function(res, input_species, source_name = NA,
                                          accepted_only = TRUE, include_synonyms = FALSE, verbose = FALSE) {
  if (is.null(res)) return(tibble(input_species = character(), source = character(), candidate_name = character(),
                                  rank = character(), status = character(), matched_name = character(),
                                  match_type = character(), outlink = character()))
  
  get_names_entries <- function(x) {
    if (is.null(x)) return(NULL)
    if ("names" %in% names(x) && length(x$names) > 0) return(x$names)
    if ("data" %in% names(x) && length(x$data) > 0) {
      all <- list()
      for (el in x$data) if (!is.null(el$names) && length(el$names) > 0) all <- c(all, el$names)
      if (length(all) > 0) return(all)
    }
    if ("results" %in% names(x) && length(x$results) > 0) {
      all <- list()
      for (el in x$results) if (!is.null(el$names) && length(el$names) > 0) all <- c(all, el$names)
      if (length(all) > 0) return(all)
    }
    recursive_find <- function(y) {
      if (is.null(y)) return(NULL)
      if (is.list(y) && "names" %in% names(y) && length(y$names) > 0) return(y$names)
      if (is.list(y)) {
        for (el in y) {
          out <- recursive_find(el)
          if (!is.null(out)) return(out)
        }
      }
      return(NULL)
    }
    recursive_find(x)
  }
  
  names_entries <- get_names_entries(res)
  if (is.null(names_entries) || length(names_entries) == 0) {
    if (verbose) message("No 'names' entries found for ", input_species)
    return(tibble(input_species = character(), source = character(), candidate_name = character(),
                  rank = character(), status = character(), matched_name = character(),
                  match_type = character(), outlink = character()))
  }
  
  rank_indicators_regex <- "(^subsp$|^subsp\\.$|^subspecies$|^subvar$|^subvar\\.$|^var$|^var\\.|^variety$|^forma$|^f\\.$|^form$|^cultivar$|^cv$|^cv\\.|^race$)"
  
  rows <- list()
  sp_tokens <- str_split(str_squish(input_species), "\\s+")[[1]]
  sp_prefix <- if (length(sp_tokens) >= 2) tolower(paste(sp_tokens[1:2], collapse = " ")) else tolower(input_species)
  
  for (entry in names_entries) {
    result_nodes <- list()
    if (!is.null(entry$results) && length(entry$results) > 0) result_nodes <- c(result_nodes, entry$results)
    if (!is.null(entry$bestResult)) result_nodes <- c(result_nodes, list(entry))
    if (!is.null(entry$names) && length(entry$names) > 0) result_nodes <- c(result_nodes, entry$names)
    if (length(result_nodes) == 0) result_nodes <- list(entry)
    
    for (node in result_nodes) {
      match_type <- NULL
      if (!is.null(node$matchType)) match_type <- node$matchType
      if (is.null(match_type) && !is.null(node$match_type)) match_type <- node$match_type
      if (is.null(match_type) && !is.null(entry$matchType)) match_type <- entry$matchType
      if (is.null(match_type) && !is.null(entry$match_type)) match_type <- entry$match_type
      
      br <- NULL
      if (!is.null(node$bestResult)) br <- node$bestResult
      if (is.null(br) && (!is.null(node$currentName) || !is.null(node$matchedName) || !is.null(node$taxonomicStatus))) br <- node
      if (is.null(br)) next
      
      status <- NULL
      if (!is.null(br$taxonomicStatus)) status <- br$taxonomicStatus
      if ((is.null(status) || !nzchar(status)) && !is.null(br$taxonomic_status)) status <- br$taxonomic_status
      status_norm <- tolower(ifelse(is.null(status), "", status))
      
      candidate <- NULL
      if (!is.null(br$currentName) && nzchar(br$currentName)) candidate <- br$currentName
      if (is.null(candidate) && !is.null(br$matchedName) && nzchar(br$matchedName)) candidate <- br$matchedName
      if (is.null(candidate) && !is.null(node$name) && nzchar(node$name)) candidate <- node$name
      if (is.null(candidate) && !is.null(entry$name) && nzchar(entry$name)) candidate <- entry$name
      if (is.null(candidate)) next
      candidate <- trimws(candidate)
      
      rank_val <- NULL
      if (!is.null(br$rank) && nzchar(br$rank)) rank_val <- br$rank
      if (is.null(rank_val) && !is.null(br$taxonRank) && nzchar(br$taxonRank)) rank_val <- br$taxonRank
      if (is.null(rank_val) && !is.null(br$rankString) && nzchar(br$rankString)) rank_val <- br$rankString
      if (is.null(rank_val) && !is.null(node$rank) && nzchar(node$rank)) rank_val <- node$rank
      if (is.null(rank_val) && !is.null(node$taxonRank) && nzchar(node$taxonRank)) rank_val <- node$taxonRank
      
      toks <- str_split(candidate, "\\s+")[[1]]
      cand_prefix <- if (length(toks) >= 2) tolower(paste(toks[1:2], collapse = " ")) else tolower(candidate)
      
      is_infraspecific <- FALSE
      if (!is.null(rank_val) && nzchar(rank_val)) {
        rlow <- tolower(rank_val)
        if (grepl(rank_indicators_regex, rlow, perl = TRUE)) is_infraspecific <- TRUE
      }
      if (!is_infraspecific && length(toks) >= 3) {
        third <- toks[3]
        if (grepl(rank_indicators_regex, tolower(third), perl = TRUE)) is_infraspecific <- TRUE
        else if (grepl("^[a-z]", third)) is_infraspecific <- TRUE
      }
      
      if (tolower(cand_prefix) == sp_prefix && !is_infraspecific) next
      if (!is_infraspecific) next
      
      if (!include_synonyms) {
        if (accepted_only && (status_norm != "accepted")) next
      }
      
      outlink <- NA_character_
      if (!is.null(node$outlink) && nzchar(node$outlink)) outlink <- node$outlink
      else if (!is.null(node$url) && nzchar(node$url)) outlink <- node$url
      else if (!is.null(entry$results) && length(entry$results) > 0) {
        for (r in entry$results) if (!is.null(r$outlink) && nzchar(r$outlink)) { outlink <- r$outlink; break }
      }
      
      matched_name <- if (!is.null(br$matchedName)) br$matchedName else NA_character_
      
      rows[[length(rows) + 1]] <- tibble(
        input_species = input_species,
        source = ifelse(is.null(source_name), NA_character_, source_name),
        candidate_name = candidate,
        rank = ifelse(is.null(rank_val), NA_character_, rank_val),
        status = ifelse(is.null(status), NA_character_, status),
        matched_name = ifelse(is.null(matched_name), NA_character_, matched_name),
        match_type = ifelse(is.null(match_type), NA_character_, match_type),
        outlink = ifelse(is.null(outlink), NA_character_, outlink)
      )
    }
  }
  
  if (length(rows) == 0) {
    return(tibble(input_species = character(), source = character(), candidate_name = character(),
                  rank = character(), status = character(), matched_name = character(),
                  match_type = character(), outlink = character()))
  }
  
  out_df <- bind_rows(rows) %>% distinct(input_species, source, candidate_name, .keep_all = TRUE)
  out_df
}
