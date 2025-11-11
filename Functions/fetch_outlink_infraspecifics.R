# Functions/fetch_outlink_infraspecifics.R
# Fetch any outlink (WFO or GRIN) and extract infraspecific-like names.
# Uses parse_wfo_json_infraspecifics() for JSON payloads and HTML heuristics for HTML pages.

library(httr)
library(rvest)
library(stringr)

.extract_taxa_from_text <- function(text) {
  if (is.null(text) || !nzchar(text)) return(character(0))
  txt <- text
  txt <- gsub('\\"', '"', txt, perl = TRUE)
  txt <- gsub('\\\\/', "/", txt, perl = TRUE)
  txt <- gsub('\\\\n', " ", txt, perl = TRUE)
  txt <- str_squish(txt)
  taxon_pattern <- paste0("([A-Z][a-z]+\\s+[a-z\\-]+(?:\\s+(?:subsp\\.?|var\\.?|subvar\\.?|forma|f\\.?|form\\b|cv\\.?|cultivar|race)\\s+[a-z\\-]+|\\s+[a-z\\-]+)(?:\\s+[A-Za-z\\.\\-\\&'()]+)?)")
  matches <- str_extract_all(txt, regex(taxon_pattern))[[1]]
  if (length(matches) == 0) return(character(0))
  matches <- unique(trimws(as.character(matches)))
  keep <- vapply(matches, function(name) {
    toks <- str_split(str_squish(name), "\\s+")[[1]]
    if (length(toks) < 3) return(FALSE)
    third <- toks[3]
    if (grepl("^(?i)(subsp\\.?|subspecies|subvar\\.?|var\\.?|variety|forma|f\\.?|form\\b|cv\\.?|cultivar|race)$", third, perl = TRUE)) return(TRUE)
    if (grepl("^[a-z]", third)) return(TRUE)
    FALSE
  }, logical(1))
  unique(matches[keep])
}

fetch_outlink_infraspecifics <- function(outlink, species_binomial = NULL, verbose = FALSE, timeout_sec = 30) {
  if (is.null(outlink) || !nzchar(outlink)) return(list(included_varieties = character(0), synonyms = character(0), page_candidates = character(0)))
  if (verbose) message("GET (Accept:text/html) ", outlink)
  req <- try(GET(outlink, user_agent("R (fetch_outlink_infraspecifics)"), accept("text/html"), timeout(timeout_sec)), silent = TRUE)
  if (inherits(req, "try-error")) return(list(included_varieties = character(0), synonyms = character(0), page_candidates = character(0)))
  ct <- tolower(ifelse(is.null(headers(req)$`content-type`), "", headers(req)$`content-type`))
  body <- content(req, as = "text", encoding = "UTF-8")
  
  if (grepl("json", ct)) {
    parsed_json <- parse_wfo_json_infraspecifics(body, species_binomial = species_binomial, verbose = verbose)
    if ((length(parsed_json$included_varieties) + length(parsed_json$synonyms) + length(parsed_json$page_candidates)) == 0) {
      pcs <- .extract_taxa_from_text(body)
      return(list(included_varieties = character(0), synonyms = character(0), page_candidates = pcs))
    }
    return(parsed_json)
  }
  
  parsed_html <- try(read_html(body), silent = TRUE)
  if (inherits(parsed_html, "try-error")) {
    pcs <- .extract_taxa_from_text(body)
    return(list(included_varieties = character(0), synonyms = character(0), page_candidates = pcs))
  }
  
  included_varieties <- character(0); synonyms <- character(0); page_candidates <- character(0)
  h_nodes <- html_nodes(parsed_html, xpath = "//h1|//h2|//h3|//h4|//h5|//h6")
  if (length(h_nodes) > 0) {
    for (hn in h_nodes) {
      title <- str_squish(html_text(hn)); title_low <- tolower(title)
      if (str_detect(title_low, "included variety")) {
        sibs <- html_nodes(hn, xpath = "following-sibling::*")
        collected <- character(0)
        for (s in sibs) {
          if (xml_name(s) %in% c("h1","h2","h3","h4","h5","h6")) break
          if (xml_name(s) %in% c("ul","ol")) collected <- c(collected, html_text(html_nodes(s, "li")))
          else {
            a_texts <- html_text(html_nodes(s, "a")); p_texts <- html_text(html_nodes(s, "p"))
            if (length(a_texts) > 0) collected <- c(collected, a_texts)
            if (length(p_texts) > 0) collected <- c(collected, p_texts)
          }
        }
        for (c in collected) {
          found <- .extract_taxa_from_text(c)
          if (length(found) > 0) included_varieties <- c(included_varieties, found) else {
            txt <- str_squish(gsub("<[^>]+>", "", c))
            if (nzchar(species_binomial) && str_detect(txt, fixed(species_binomial, ignore_case = TRUE))) included_varieties <- c(included_varieties, txt)
          }
        }
      } else if (str_detect(title_low, "synonym")) {
        sibs <- html_nodes(hn, xpath = "following-sibling::*")
        collected <- character(0)
        for (s in sibs) {
          if (xml_name(s) %in% c("h1","h2","h3","h4","h5","h6")) break
          if (xml_name(s) %in% c("ul","ol")) collected <- c(collected, html_text(html_nodes(s, "li")))
          else {
            a_texts <- html_text(html_nodes(s, "a")); p_texts <- html_text(html_nodes(s, "p"))
            if (length(a_texts) > 0) collected <- c(collected, a_texts)
            if (length(p_texts) > 0) collected <- c(collected, p_texts)
          }
        }
        for (c in collected) {
          found <- .extract_taxa_from_text(c)
          if (length(found) > 0) synonyms <- c(synonyms, found) else {
            txt <- str_squish(gsub("<[^>]+>", "", c))
            if (nzchar(species_binomial) && str_detect(txt, fixed(species_binomial, ignore_case = TRUE))) synonyms <- c(synonyms, txt)
          }
        }
      }
    }
  }
  
  page_candidates <- unique(c(page_candidates, .extract_taxa_from_text(html_text(parsed_html))))
  list(included_varieties = unique(included_varieties), synonyms = unique(synonyms), page_candidates = unique(page_candidates))
}
