# Infraspecifics pipeline — README

What it does
- Queries the Global Names Verifier (WFO + GRIN), follows verifier outlinks, extracts infraspecific names from WFO/HTML pages, cleans and re-verifies candidates, and produces a slim long table of verified infraspecific names with source attribution.

Required files
- scripts/compile_reverify_infraspecifics_full_with_write.R  — main pipeline
- Functions/Query_taxa_resolver.R
- Functions/Process_taxa_resolver_results.R
- Functions/parse_wfo_json_infraspecifics.R
- Functions/fetch_outlink_infraspecifics.R
- Tools/clean_outlink_tokens.R

Dependencies
- R packages: httr, jsonlite, stringr, dplyr, tidyr, readr, writexl, tibble
- Install in R: install.packages(c("httr","jsonlite","stringr","dplyr","tidyr","readr","writexl","tibble"))

Configuration (quick)
- Edit species_list inside scripts/compile_reverify_infraspecifics_full_with_write.R to set target species.
- Adjust query_delay and outlink_delay in the script to be polite to APIs.

Output
- infraspecifics_verified.xlsx — one row per infraspecific with columns:
  - input_species — queried species name
  - infraspecific_name — each verified infraspecific name (one per row)
  - source — provenance: "WFO", "GRIN", "WFO, GRIN" (both), "OUTLINK", or "NONE"

Source & cache
- Per-name source checks use the verifier and are cached to: verifier_provenance_cache.rds
- To force fresh provenance checks delete that file before rerunning.

Quick hints
- To inspect raw parsing for a species (interactive): examine the in-memory object all_results_per_species after running the script (contains raw matches and parsed outlink payloads).
- To tune filtering: edit noise_words, infra_syntax_ok, or accept_infra_regex in the script.
- For many species, increase query_delay and outlink_delay to avoid rate limits and consider persisting intermediate results.

Troubleshooting
- If the final Excel is empty, check that helper modules are present and that API responses are reachable from your environment.
- If noisy tokens persist, add problematic words to the noise_words list in the script.
