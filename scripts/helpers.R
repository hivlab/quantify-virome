
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(magrittr)
library(xml2)
library(httr)
library(rvest)
library(readr)

# Split concatenated BLAST+ xml
split_hits <- function(xml) {
  n <- diff(c(which(str_detect(xml, "<\\?xml")), length(xml) + 1))
  split(xml, rep(1:length(n), n))
}

# Parse BLAST+ xml
parse_hits <- function(xml) {
  program <- xml_find_first(xml, "/BlastOutput//BlastOutput_program") %>% xml_text()
  db <- xml_find_first(xml, "/BlastOutput//BlastOutput_db") %>% xml_text() %>% basename()
  iter <- xml_find_all(xml, "/BlastOutput//Iteration")
  query <- xml_find_first(iter, "Iteration_query-def") %>% xml_text()
  hits <- xml_find_all(iter, "Iteration_hits")

  children <- xml_children(hits) %>%
    map(xml_children)
  gi <- map(children, 2) %>%
    map(xml_text) %>%
    unlist() %>%
    str_replace("gi\\|(\\d+)\\|.*", "\\1")
  Hit_accession <- map(children, 4) %>% map(xml_text) %>% unlist()
  Hit_def <- map(children, 3) %>% map(xml_text) %>% unlist()

  hsps <- map(children, 6) %>% map(xml_children) %>% map(1) %>% map(xml_children)
  hsps_txt <- map(hsps, xml_text)
  hsps_nm <- map(hsps, xml_name)
  hsps <- map2(hsps_nm, hsps_txt, ~tibble(key = .x, value = .y)) %>%
    map(spread, key, value) %>%
    bind_rows() %>%
    mutate_at(vars(`Hsp_align-len`, `Hsp_bit-score`, `Hsp_evalue`, `Hsp_gaps`, `Hsp_hit-frame`,
                   `Hsp_hit-from`, `Hsp_hit-to`, Hsp_identity, Hsp_num, Hsp_positive,
                   `Hsp_query-frame`, `Hsp_query-from`, `Hsp_query-to`, Hsp_score), parse_number)

  data_frame(program, db, query, gi, Hit_accession, Hit_def) %>% bind_cols(hsps)
}

query_taxid <- function(gi) {
  res <- GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
             query = list(db = "nucleotide", id = gi, rettype = "fasta", retmode = "xml"))
  cont <- content(res, as = "parsed", encoding = "UTF-8")
  xml_node(cont, xpath = "//TSeq_taxid") %>% xml_text()
}

#' @param blast_xml Path to BLAST+ XML result file, a character string.
#' @param ... further path(s) to BLAST+ XML result file(s), a character string.
parse_blast_xml <- function(blast_xml, ...) {

  message("Open and fix BLAST+ xml strings")
  xml_str <- data_frame(path = c(blast_xml, ...)) %>%
    mutate(xml = map(path, read_lines)) %>%
    mutate(xml = map(xml, split_hits)) %>%
    unnest() %>%
    mutate(xml = map(xml, str_c, collapse = ""))

  message("Import fixed xml")
  xml_imported <- mutate(xml_str, xml = map(xml, read_xml))

  message("Extract hits from xml")
  xml_parsed <- mutate(xml_imported, hits = map(xml, parse_hits))

  message("Extract blast results from list")
  blast_results <- dplyr::select(xml_parsed, hits) %>% unnest()
  class(blast_results) <- append(class(blast_results), "blast_xml_tab")
  return(blast_results)
}

gi2taxid <- function(tab, taxdb) {
  UseMethod("gi2taxid", tab)
}

#' @param tab data_frame with parsed blast results
#' @param taxdb sqlite database with gi_taxid_nucl and gi_taxid_prot tables
gi2taxid.blast_xml_tab <- function(tab, taxdb) {

  mapped_gis <- tab$gi

  message("Connect to database")
  db <- src_sqlite(taxdb, create = TRUE)

  message("Collect tax_ids from tables")
  gi_nuc <- tbl(db, "gi_taxid_nucl") %>%
    filter(gi %in% mapped_gis) %>%
    collect()
  gi_prot <- tbl(db, "gi_taxid_prot") %>%
    filter(gi %in% mapped_gis) %>%
    collect()

  message("Bind rows")
  gi_tab <- bind_rows(gi_nuc, gi_prot) %>%
    mutate_at(vars(gi), as.character)

  message("Join taxonomy to blast results by gi")
  known <- left_join(tab, gi_tab)

  message("Fill in few missing tax_ids by quering remote ncbi database")
  with_taxid <- filter(known, !is.na(tax_id))
  no_taxid <- filter(known, is.na(tax_id))

  query <- mutate(no_taxid, tax_id = map_chr(gi, query_taxid)) %>%
    dplyr::select(gi, tax_id) %>%
    mutate_at("tax_id", as.integer)

  fixed_taxid <- full_join(dplyr::select(no_taxid, -tax_id), query)
  blast_results_taxids <- bind_rows(with_taxid, fixed_taxid)
  class(blast_results_taxids) <- append(class(blast_results_taxids), "blast_results_taxids")
  return(blast_results_taxids)
}

filter_division <- function(tab, nodes, div_id, div, not_div) {
  UseMethod("filter_division", tab)
}

#' @param tab blast results tab with tax_ids, a data_frame.
#' @param nodes path to nodes.csv file, a character string.
#' @param div_id division id of interest, integer. Defaults to 3, viruses.
#' @param div path ot output.csv file for records belonging to div_id, a character string.
#' @param not_div path ot output.csv file for records NOT belongigng to div_id, a character string.
filter_division.blast_results_taxids <- function(tab, nodes, div_id, div, not_div) {

  message("Import gi tax_id table and taxonomy nodes table")
  nodes <- read_csv(nodes)

  message("Merge names by tax_id to known sequences")
  mapped_tab <- left_join(tab, nodes)

  division <- filter(mapped_tab, near(division_id, div_id))
  not_division <- filter(mapped_tab, !near(division_id, div_id))

  write_csv(division, div)
  write_csv(not_division, not_div)
}
