# End-to-end pipeline
# Extract infraspecific taxon names from World Flora Online and GRIN-Global
# query global names verifier, parse outlinks, clean, re-verify, expand to long table
# Output: single Excel with columns (input_species, infraspecific_name, source)
#
# Notes:
# - Run after sourcing Functions/ Tools/ and helpers in the project.
# - Produces one Excel file: infraspecifics_verified.xlsx


# ---- dependencies ----
required_pkgs <- c("httr","jsonlite","stringr","dplyr","tidyr","rvest","xml2","readr","writexl","tibble")
for (p in required_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(httr)
library(jsonlite)
library(stringr)
library(dplyr)
library(tidyr)
library(rvest)
library(xml2)
library(readr)
library(writexl)
library(tibble)

# ---- source helper modules (adjust paths if needed) ----
source('Functions/Query_taxa_resolver.R')
source('Functions/Process_taxa_resolver_results.R')
source('Functions/parse_wfo_json_infraspecifics.R')
source('Functions/fetch_outlink_infraspecifics.R')
source('Tools/clean_outlink_tokens.R')

# ---- configuration with example species list ----
species_list <- c("Capsicum annuum", "Cynara cardunculus", "Brassica oleracea")
sources <- list(WFO = '196', GRIN = '6')
query_delay <- 0.12
outlink_delay <- 0.18

# ---- helper functions ----
# build a small set of query variants for each species
expand_queries_for_species <- function(species) {
  toks <- str_split(str_squish(species), "\\s+")[[1]]
  if (length(toks) < 2) return(character(0))
  genus <- toks[1]
  rank_prefixes <- c("var.", "variety", "subsp.", "subspecies", "f.", "forma", "subvar.", "subvariety", "cv.", "cultivar")
  unique(c(genus, species, paste(species, rank_prefixes), paste0(species, " ")))
}

# canonicalize an infraspecific to a short key (genus species rank epithet) or NA
canonicalize_infra <- function(name) {
  if (is.na(name) || !nzchar(name)) return(NA_character_)
  nm <- str_squish(as.character(name))
  toks_raw <- strsplit(nm, "\\s+")[[1]]
  toks_clean <- tolower(gsub("[[:punct:]]", "", toks_raw))
  if (length(toks_clean) < 3) return(NA_character_)
  rank_tokens <- c("subsp","subspecies","subvar","var","variety","forma","form","f","cultivar","cv","race")
  rank_idx <- which(toks_clean %in% rank_tokens)
  if (length(rank_idx) >= 1) {
    idx <- rank_idx[1]
    if (length(toks_clean) >= idx + 1) {
      return(paste(toks_clean[1], toks_clean[2], toks_clean[idx], toks_clean[idx + 1], sep = " "))
    }
  }
  third <- toks_raw[3]
  if (grepl("^[a-z]", third)) return(paste(toks_clean[1], toks_clean[2], toks_clean[3], sep = " "))
  NA_character_
}

# local wrapper for verifier queries (WFO+GRIN by default)
query_taxa_resolver_local <- function(taxa, sources_ids = c('196','6'), timeout_sec = 30) {
  if (!is.character(taxa) || length(taxa) != 1) return(NULL)
  taxa_encoded <- URLencode(taxa, reserved = TRUE)
  URL <- paste0('https://verifier.globalnames.org/api/v1/verifications/', taxa_encoded,
                '?data_sources=', paste(sources_ids, collapse = "|"),
                '&all_matches=true&capitalize=true&species_group=true',
                '&fuzzy_uninomial=false&stats=false&main_taxon_threshold=0.8')
  res <- try(GET(URL, user_agent("R (reverify)"), timeout(timeout_sec)), silent = TRUE)
  if (inherits(res, "try-error") || status_code(res) != 200) return(NULL)
  txt <- content(res, as = "text", encoding = "UTF-8")
  parsed <- try(fromJSON(txt, simplifyVector = FALSE), silent = TRUE)
  if (inherits(parsed, "try-error")) return(NULL)
  parsed
}

# quick predicate: does a candidate look like an infraspecific for given species?
infra_syntax_ok <- function(name, species_binomial) {
  if (is.null(name) || !nzchar(name)) return(FALSE)
  genus <- str_split(str_squish(species_binomial), "\\s+")[[1]][1]
  pattern <- paste0("^\\s*(", genus, "|", tolower(genus), ")\\s+[a-z\\-]+(?:\\s+(?:subsp\\.?|subspecies|subvar\\.?|var\\.?|variety|forma|f\\.?|form\\b|cv\\.?|cultivar|race)\\s+[a-z\\-]+|\\s+[a-z\\-]+)(?:\\s+.*)?$")
  if (!str_detect(name, regex(pattern))) return(FALSE)
  toks <- str_split(str_squish(name), "\\s+")[[1]]
  if (length(toks) < 3) return(FALSE)
  third <- toks[3]
  noise_words <- c("and","or","of","in","the","using","application","implications","exploitable","genotypes","estimated","estimate","assessment","analysis","based","revealed")
  if (tolower(third) %in% noise_words) return(FALSE)
  if (!grepl("^[a-z\\-]+$", third)) {
    if (!grepl("^(subsp\\.?|subspecies|subvar\\.?|var\\.?|variety|forma|f\\.?|form\\b|cv\\.?|cultivar|race)$", third, perl = TRUE)) {
      return(FALSE)
    }
  }
  TRUE
}

# ---- PHASE 1: verifier expansion and outlink collection ----
all_results_per_species <- list()
raw_summary_rows <- list()

for (sp in species_list) {
  queries <- expand_queries_for_species(sp)
  combined_matches_list <- list()
  raw_responses <- list()
  
  for (q in queries) {
    for (src_label in names(sources)) {
      Sys.sleep(query_delay)
      src_id <- sources[[src_label]]
      res <- query_taxa_resolver(q, sources = c(src_id), verbose = FALSE)
      raw_responses[[paste(q, src_label, sep = "|||")]] <- res
      df_matches <- tryCatch({
        extract_infraspecifics_single(res, input_species = sp, source_name = src_label,
                                      accepted_only = FALSE, include_synonyms = TRUE, verbose = FALSE)
      }, error = function(e) tibble())
      if (nrow(df_matches) > 0) combined_matches_list[[length(combined_matches_list) + 1]] <- df_matches
    }
  }
  
  combined_matches <- if (length(combined_matches_list) == 0) {
    tibble(input_species = character(), source = character(), candidate_name = character(),
           rank = character(), status = character(), matched_name = character(), match_type = character(), outlink = character())
  } else bind_rows(combined_matches_list) %>% distinct(input_species, source, candidate_name, .keep_all = TRUE)
  
  # extract outlinks
  outlinks <- unique(na.omit(unlist(lapply(raw_responses, function(x) {
    if (is.null(x) || is.null(x$names)) return(NULL)
    ols <- c()
    for (n in x$names) {
      if (!is.null(n$results)) {
        for (r in n$results) if (!is.null(r$outlink) && nzchar(r$outlink)) ols <- c(ols, r$outlink)
      }
      if (!is.null(n$outlink) && nzchar(n$outlink)) ols <- c(ols, n$outlink)
    }
    ols
  }))))
  
  # follow outlinks and parse
  parsed_from_outlinks <- list()
  if (length(outlinks) > 0) {
    for (ol in outlinks) {
      Sys.sleep(outlink_delay)
      parsed <- tryCatch(fetch_outlink_infraspecifics(ol, species_binomial = sp, verbose = FALSE), error = function(e) NULL)
      if (!is.null(parsed)) parsed_from_outlinks[[ol]] <- parsed
    }
  }
  
  # collect parsed names
  wfo_included <- unique(unlist(lapply(parsed_from_outlinks, function(x) if (!is.null(x$included_varieties)) x$included_varieties else character(0))))
  wfo_syns     <- unique(unlist(lapply(parsed_from_outlinks, function(x) if (!is.null(x$synonyms)) x$synonyms else character(0))))
  page_candidates <- unique(unlist(lapply(parsed_from_outlinks, function(x) if (!is.null(x$page_candidates)) x$page_candidates else character(0))))
  
  # append parsed names
  appended <- combined_matches
  if (length(wfo_included) > 0) {
    appended <- bind_rows(appended,
                          tibble(input_species = sp, source = "WFO_OUTLINK", candidate_name = wfo_included,
                                 rank = NA_character_, status = NA_character_, matched_name = NA_character_, match_type = NA_character_, outlink = NA_character_))
  }
  if (length(wfo_syns) > 0) {
    appended <- bind_rows(appended,
                          tibble(input_species = sp, source = "WFO_OUTLINK", candidate_name = wfo_syns,
                                 rank = NA_character_, status = NA_character_, matched_name = NA_character_, match_type = NA_character_, outlink = NA_character_))
  }
  if (length(page_candidates) > 0) {
    appended <- bind_rows(appended,
                          tibble(input_species = sp, source = "OUTLINK_PAGE", candidate_name = page_candidates,
                                 rank = NA_character_, status = NA_character_, matched_name = NA_character_, match_type = NA_character_, outlink = NA_character_))
  }
  
  appended <- appended %>% distinct(input_species, source, candidate_name, .keep_all = TRUE)
  
  # aggregated raw row
  wfo_cands <- appended %>% filter(source %in% c("WFO","WFO_OUTLINK")) %>% pull(candidate_name) %>% unique()
  grin_cands <- appended %>% filter(source %in% c("GRIN","GRIN_OUTLINK")) %>% pull(candidate_name) %>% unique()
  outlink_page_cands <- appended %>% filter(source == "OUTLINK_PAGE") %>% pull(candidate_name) %>% unique()
  
  raw_row <- data.frame(
    input_species = sp,
    infraspecific_list_GRIN = if (length(grin_cands) == 0) NA_character_ else paste(sort(unique(grin_cands)), collapse = "; "),
    n_infraspecific_GRIN = length(unique(grin_cands)),
    infraspecific_list_WFO = if (length(wfo_cands) == 0) NA_character_ else paste(sort(unique(wfo_cands)), collapse = "; "),
    n_infraspecific_WFO = length(unique(wfo_cands)),
    infraspecific_list_WFO_included_variety = if (length(wfo_included) == 0) NA_character_ else paste(sort(unique(wfo_included)), collapse = "; "),
    n_infraspecific_WFO_included_variety = length(unique(wfo_included)),
    infraspecific_list_WFO_synonyms = if (length(wfo_syns) == 0) NA_character_ else paste(sort(unique(wfo_syns)), collapse = "; "),
    n_infraspecific_WFO_synonyms = length(unique(wfo_syns)),
    infraspecific_list_OUTLINK_PAGE = if (length(outlink_page_cands) == 0) NA_character_ else paste(sort(unique(outlink_page_cands)), collapse = "; "),
    n_infraspecific_OUTLINK_PAGE = length(unique(outlink_page_cands)),
    stringsAsFactors = FALSE
  )
  
  raw_summary_rows[[length(raw_summary_rows) + 1]] <- raw_row
  all_results_per_species[[sp]] <- list(raw_matches = appended, outlink_parsed = parsed_from_outlinks)
}

raw_summary_df <- bind_rows(raw_summary_rows)

# ---- PHASE 2: clean OUTLINK_PAGE tokens (genus-aware) ----
cleaned_df <- clean_outlinks_in_df(raw_summary_df)

# ---- PHASE 3: re-verify extracted combined candidates (stricter) ----
collect_candidates_for_reverify <- function(row) {
  sp <- row$input_species
  cells <- c(row$infraspecific_list_GRIN, row$infraspecific_list_WFO,
             row$infraspecific_list_WFO_included_variety, row$infraspecific_list_WFO_synonyms,
             row$infraspecific_list_OUTLINK_PAGE_clean, row$infraspecific_list_combined)
  parts <- unique(unlist(lapply(cells, function(cell) {
    if (is.na(cell) || !nzchar(cell)) return(character(0))
    unlist(str_split(as.character(cell), "\\s*;\\s*"))
  })))
  genus <- str_split(sp, "\\s+")[[1]][1]
  parts <- unique(parts[grepl(paste0("^", genus, "\\b"), parts)])
  parts
}

reverified_summary_rows <- list()

for (i in seq_len(nrow(cleaned_df))) {
  row <- cleaned_df[i, ]
  sp <- row$input_species
  cands <- collect_candidates_for_reverify(row)
  verified_names <- character(0)
  if (length(cands) > 0) {
    for (cand in cands) {
      if (!infra_syntax_ok(cand, sp)) next
      Sys.sleep(query_delay)
      res <- query_taxa_resolver_local(cand, sources_ids = c(sources$WFO, sources$GRIN))
      if (is.null(res)) next
      names_entries <- NULL
      if (!is.null(res$names)) names_entries <- res$names
      else if (!is.null(res$data)) names_entries <- unlist(lapply(res$data, function(d) d$names), recursive = FALSE)
      else if (!is.null(res$results)) names_entries <- unlist(lapply(res$results, function(r) r$names), recursive = FALSE)
      if (is.null(names_entries) || length(names_entries) == 0) next
      for (ne in names_entries) {
        br <- if (!is.null(ne$bestResult)) ne$bestResult else ne
        match_type <- NULL
        if (!is.null(br$matchType)) match_type <- br$matchType
        if (is.null(match_type) && !is.null(br$match_type)) match_type <- br$match_type
        if (is.character(match_type) && tolower(as.character(match_type)) == "nomatch") next
        cand_name <- if (!is.null(br$currentName)) br$currentName else if (!is.null(br$matchedName)) br$matchedName else if (!is.null(ne$name)) ne$name else NA_character_
        if (is.na(cand_name)) next
        if (!infra_syntax_ok(cand_name, sp)) next
        if (is.na(canonicalize_infra(cand_name))) next
        verified_names <- c(verified_names, cand_name)
      }
    }
  }
  verified_names <- unique(verified_names)
  reverified_summary_rows[[length(reverified_summary_rows) + 1]] <- data.frame(
    input_species = sp,
    infraspecific_list_verified = if (length(verified_names) == 0) NA_character_ else paste(sort(unique(verified_names)), collapse = "; "),
    n_infraspecific_verified = length(unique(verified_names)),
    stringsAsFactors = FALSE
  )
}

reverified_df <- bind_rows(reverified_summary_rows)
final_output_df <- cleaned_df %>% left_join(reverified_df, by = "input_species")

# ---- PHASE 4: final cleaning of verified names ----
noise_words <- c(
  "and", "or", "of", "in", "the", "using", "application", "implications", "exploitable",
  "genotypes", "genebank", "n", "species with n", "estimated", "estimate", "assessment", "analysis", "based", "revealed",
  "using", "names", "accepted", "click", "please", "exportable", "if", "load", "requests",
  "ensure", "check", "snapshot", "atlas", "catalogue", "button", "view", "find", "explore", "choose"
)
noise_rx <- paste0("\\b(", paste(noise_words, collapse = "|"), ")\\b")

accept_infra_regex <- function(species_binomial) {
  gen <- str_split(str_squish(species_binomial), "\\s+")[[1]][1]
  paste0("^\\s*(", gen, "|", tolower(gen), ")\\s+[a-z\\-]+(?:\\s+(?:subsp\\.?|subspecies|subvar\\.?|var\\.?|variety|forma|f\\.?|form\\b|cv\\.?|cultivar|race)\\s+[a-z\\-]+|\\s+[a-z\\-]+)(?:\\s+.*)?$")
}

clean_verified_cell <- function(verified_cell, species_binomial) {
  if (is.na(verified_cell) || !nzchar(verified_cell)) return(NA_character_)
  parts <- unlist(str_split(as.character(verified_cell), "\\s*;\\s*"))
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  if (length(parts) == 0) return(NA_character_)
  rx <- accept_infra_regex(species_binomial)
  keep <- vapply(parts, function(p) {
    toks <- str_split(str_squish(p), "\\s+")[[1]]
    if (length(toks) > 0 && tolower(toks[1]) %in% noise_words) return(FALSE)
    if (str_detect(tolower(p), noise_rx)) {
      if (length(toks) <= 3) return(FALSE)
      if (any(tolower(toks[1:min(3, length(toks))]) %in% noise_words)) return(FALSE)
    }
    if (!str_detect(p, regex(rx))) return(FALSE)
    if (is.na(canonicalize_infra(p))) return(FALSE)
    TRUE
  }, logical(1))
  kept <- unique(parts[keep])
  if (length(kept) == 0) return(NA_character_)
  paste(kept, collapse = "; ")
}

final_df <- final_output_df %>% rowwise() %>%
  mutate(
    infraspecific_list_verified_clean = clean_verified_cell(infraspecific_list_verified, input_species),
    n_infraspecific_verified_clean = ifelse(is.na(infraspecific_list_verified_clean), 0L, lengths(strsplit(infraspecific_list_verified_clean, "\\s*;\\s*")))
  ) %>% ungroup()

# ---- PHASE 5: expand to one-row-per-infraspecific and determine source ----
prov_cache_file <- "verifier_provenance_cache.rds"
if (file.exists(prov_cache_file)) {
  cached_list <- readRDS(prov_cache_file)
  if (is.list(cached_list)) {
    prov_cache <- list2env(cached_list, parent = emptyenv())
  } else {
    prov_cache <- new.env(parent = emptyenv())
  }
} else prov_cache <- new.env(parent = emptyenv())

query_source_for_name <- function(name, source_id, timeout_sec = 30) {
  cache_key <- paste0("src:", source_id, "|", name)
  if (exists(cache_key, envir = prov_cache, inherits = FALSE)) return(get(cache_key, envir = prov_cache))
  res <- tryCatch(query_taxa_resolver(name, sources = c(as.character(source_id)), verbose = FALSE, timeout_sec = timeout_sec),
                  error = function(e) NULL)
  out <- list(ok = FALSE, matched = NA_character_, matchType = NA_character_, status = NA_character_)
  if (!is.null(res)) {
    names_entries <- NULL
    if (!is.null(res$names)) names_entries <- res$names
    else if (!is.null(res$data)) names_entries <- unlist(lapply(res$data, function(d) d$names), recursive = FALSE)
    else if (!is.null(res$results)) names_entries <- unlist(lapply(res$results, function(r) r$names), recursive = FALSE)
    if (!is.null(names_entries) && length(names_entries) > 0) {
      for (ne in names_entries) {
        br <- if (!is.null(ne$bestResult)) ne$bestResult else ne
        mt <- NULL
        if (!is.null(br$matchType)) mt <- br$matchType
        if (is.null(mt) && !is.null(br$match_type)) mt <- br$match_type
        if (is.character(mt) && tolower(as.character(mt)) == "nomatch") next
        cand_name <- if (!is.null(br$currentName)) br$currentName else if (!is.null(br$matchedName)) br$matchedName else if (!is.null(ne$name)) ne$name else NA_character_
        status <- if (!is.null(br$taxonomicStatus)) br$taxonomicStatus else if (!is.null(br$taxonomic_status)) br$taxonomic_status else NA_character_
        if (!is.na(cand_name)) {
          out$ok <- TRUE
          out$matched <- cand_name
          out$matchType <- mt
          out$status <- status
          break
        }
      }
    }
  }
  assign(cache_key, out, envir = prov_cache)
  try(saveRDS(as.list(prov_cache), prov_cache_file), silent = TRUE)
  out
}

determine_source <- function(candidate) {
  key_master <- paste0("prov|", candidate)
  if (exists(key_master, envir = prov_cache, inherits = FALSE)) return(get(key_master, envir = prov_cache))
  wfo <- query_source_for_name(candidate, "196")
  grin <- query_source_for_name(candidate, "6")
  src <- "NONE"
  if (isTRUE(wfo$ok) && isTRUE(grin$ok)) src <- "BOTH"
  else if (isTRUE(wfo$ok)) src <- "WFO"
  else if (isTRUE(grin$ok)) src <- "GRIN"
  else src <- "OUTLINK"
  res <- list(source = src, wfo = wfo, grin = grin)
  assign(key_master, res, envir = prov_cache)
  try(saveRDS(as.list(prov_cache), prov_cache_file), silent = TRUE)
  res
}

# Build long-form table (one row per infraspecific)
long_rows <- list()
verified_col_prefer <- if ("infraspecific_list_verified_clean" %in% names(final_df)) "infraspecific_list_verified_clean" else "infraspecific_list_verified"

for (i in seq_len(nrow(final_df))) {
  sp <- as.character(final_df$input_species[i])
  verified_cell <- final_df[[verified_col_prefer]][i]
  candidates <- character(0)
  if (!is.na(verified_cell) && nzchar(as.character(verified_cell))) {
    candidates <- unlist(str_split(as.character(verified_cell), "\\s*;\\s*"))
    candidates <- trimws(candidates)
    candidates <- candidates[candidates != ""]
    candidates <- unique(candidates)
  }
  if (length(candidates) == 0) next
  for (nm in candidates) {
    prov <- determine_source(nm)
    matched <- NA_character_; status <- NA_character_
    if (!is.null(prov$wfo) && isTRUE(prov$wfo$ok)) {
      matched <- prov$wfo$matched; status <- prov$wfo$status
    } else if (!is.null(prov$grin) && isTRUE(prov$grin$ok)) {
      matched <- prov$grin$matched; status <- prov$grin$status
    }
    long_rows[[length(long_rows) + 1]] <- data.frame(
      input_species = sp,
      infraspecific_name = nm,
      source = prov$source,
      verifier_matched_name = ifelse(is.null(matched), NA_character_, matched),
      verifier_status = ifelse(is.null(status), NA_character_, status),
      stringsAsFactors = FALSE
    )
    Sys.sleep(0.06)
  }
}

if (length(long_rows) > 0) long_df <- bind_rows(long_rows) %>% distinct(input_species, infraspecific_name, .keep_all = TRUE)
else long_df <- tibble(input_species = character(), infraspecific_name = character(), source = character(),
                       verifier_matched_name = character(), verifier_status = character())

# ---- SLIM the long_df: keep only input_species, infraspecific_name, source
#      and replace "BOTH" with "WFO, GRIN" ----
long_slim <- long_df %>%
  mutate(source = ifelse(toupper(trimws(source)) == "BOTH", "WFO, GRIN", source)) %>%
  select(input_species, infraspecific_name, source) %>%
  distinct(input_species, infraspecific_name, source)

# Write ONLY this Excel file as the final output
final_excel <- "infraspecifics_verified.xlsx"
writexl::write_xlsx(long_slim, final_excel)

##### end script ######
