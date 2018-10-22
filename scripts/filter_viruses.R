
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

#' @param ... Path(s) to BLAST+ XML result file(s), a character string.
parse_blast_xml <- function(...) {

  message("Open and fix BLAST+ xml strings")
  xml_lines <- data_frame(path = c(...)) %>%
    mutate(xml = map(path, read_lines)) %>%
    filter(!map_lgl(xml, identical, character(0)))

  if(nrow(xml_lines) == 0) stop("Empty file(s). No BLAST hits to parse!")

  xml_str <- mutate(xml_lines, xml = map(xml, split_hits)) %>%
    unnest() %>%
    mutate(xml = map(xml, str_c, collapse = ""))

  message("Import fixed xml")
  xml_imported <- mutate(xml_str, xml = map(xml, read_xml))

  message("Extract hits from xml")
  xml_parsed <- mutate(xml_imported, hits = map(xml, parse_hits))

  message("Extract blast results from list")
  blast_results <- dplyr::select(xml_parsed, hits) %>% unnest()
  class(blast_results) <- append(class(blast_results), "blast_tab")
  return(blast_results)
}

# create empty blast_xml_tab when no blast hits in sample
empty_blast_xml_tab <- function() {
  tab <- data_frame(program = character(0), db = character(0), query = character(0), gi = character(0),
                    Hit_accession = character(0), Hit_def = character(0), `Hsp_align-len` = numeric(0),
                    `Hsp_bit-scorev` = numeric(0),
                    Hsp_evalue = numeric(0), Hsp_gaps = numeric(0), `Hsp_hit-frame` = numeric(0), `Hsp_hit-from` = numeric(0),
                    `Hsp_hit-to` = numeric(0), Hsp_hseq = character(0), Hsp_identity = numeric(0), Hsp_midline = character(0),
                    Hsp_num = numeric(0), Hsp_positive = numeric(0), Hsp_qseq = character(0), `Hsp_query-frame` = numeric(0),
                    `Hsp_query-from` = numeric(0), `Hsp_query-to` = numeric(0), Hsp_score = numeric(0))
  class(tab) <- append(class(tab), "blast_tab")
  return(tab)
}

# put safe wrapping around parse_blast function and return empty blast_xml_tab when no blast hits in sample
parse_blast_xml_safe <- safely(parse_blast_xml, otherwise = empty_blast_xml_tab())


parse_blast_tsv <- function(...) {

  message("Importing BLAST+ tsv")
  blast_results <- data_frame(path = c(...)) %>%
    mutate(tsv = map(path, read_tsv)) %>%
    unnest()

  message("Munging metadata")
  blast_results <- rename(blast_results, query = qseqid, gi = sgi) %>%
    mutate_at(vars(gi), as.character) %>%
    mutate(path = basename(path)) %>%
    select(path, query, gi, everything())
  class(blast_results) <- append(class(blast_results), "blast_tab")
  return(blast_results)
}

empty_blast_tsv_tab <- function() {
  tab <- data_frame(path = character(0), query = character(0), gi = character(0))
  class(tab) <- append(class(tab), "blast_tab")
  return(tab)
}

parse_blast_tsv_safe <- safely(parse_blast_tsv, otherwise = empty_blast_tsv_tab())

gi2taxid <- function(tab, taxdb) {
  UseMethod("gi2taxid", tab)
}

#' @param tab data_frame with parsed blast results
#' @param taxdb sqlite database with gi_taxid_nucl and gi_taxid_prot tables
gi2taxid.blast_tab <- function(tab, taxdb) {

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

#' @param blast_results_taxids blast results tab with tax_ids, a data_frame.
#' @param nodes path to nodes.csv file, a character string.

#' @param div path ot output.csv file for records belonging to div_id, a character string.
#' @param not_div path ot output.csv file for records NOT belongigng to div_id, a character string.
filter_division.blast_results_taxids <- function(blast_results_taxids, nodes, div_id, div, not_div) {

  message("Import gi tax_id table and taxonomy nodes table")
  nodes <- read_csv(nodes)

  message("Merge names by tax_id to known sequences")
  mapped_tab <- left_join(blast_results_taxids, nodes)

  division <- filter(mapped_tab, division_id == div_id)
  not_division <- filter(mapped_tab, division_id != div_id)

  write_csv(division, div)
  write_csv(not_division, not_div)
}

#' Parse and add taxonomy to blast+ results
#' @param ... paths(s) to blast results xml- or tsv file, a chracter string.
#' @param div_id division id of interest, integer. Defaults to 3, viruses.
#' @param taxdb path to taxonomy database, sqlite database, a character string.
#' @param nodes path to taxonomy nodes.csv file, a chracter string.
#' @param phages path to output file for phage hits, a character string.
#' @param viruses path to output file for virus hits, character string
#' @div_id taxonomy divition id, integer.
#'
blast_taxonomy <- function(..., taxdb, nodes, phages, viruses, div_id = 3) {

  message("Parse blast results\n")
  dots <- c(...)

  if (any(grepl("xml$", dots))) {

    blast_results <- grep("xml$", dots, value = TRUE)
    tab <- parse_blast_xml_safe(blast_results)

  } else if (any(grepl("tsv$", dots))) {

    blast_results <- grep("tsv$", dots, value = TRUE)
    tab <- parse_blast_tsv_safe(blast_results)

  } else {

    stop("No files with blast results!")

  }

  # if it's not empty result then stop error
  if (!is.null(tab$error) && !str_detect(tab$error$message, "No BLAST hits to parse")) stop(tab$error)

  message("Map tax_ids to gis\n")
  tab <- gi2taxid(tab$result, taxdb)

  message("Filter phages (division_id == 3) and save results to csv files\n")
  filter_division(blast_results_taxids = tab, nodes = nodes, div = phages, not_div = viruses, div_id = div_id)
}

do.call(blast_taxonomy, c(snakemake@input, snakemake@output))

