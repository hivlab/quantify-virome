
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(magrittr)
library(xml2)
library(httr)
library(rvest)
library(readr)

query_taxid <- function(gi) {
  res <- httr::GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
             query = list(db = "nucleotide", id = gi, rettype = "fasta", retmode = "xml"))
  cont <- httr::content(res, as = "parsed", encoding = "UTF-8")
  rvest::xml_node(cont, xpath = "//TSeq_taxid") %>% xml2::xml_text()
}

parse_blast_tsv <- function(...) {

  message("Importing BLAST+ tabular output")
  blast_results <- dplyr::data_frame(path = c(...)) %>%
    dplyr::dplyr::mutate(tsv = purrr::map(path, readr::read_tsv)) %>%
    tidyr::unnest()

  message("Munging metadata")
  blast_results <- dplyr::rename(blast_results, query = qseqid, gi = sgi) %>%
    dplyr::mutate_at(dplyr::vars(gi), as.character) %>%
    dplyr::mutate(path = basename(path)) %>%
    dplyr::select(path, query, gi, dplyr::everything())
  class(blast_results) <- append(class(blast_results), "blast_tabular")
  return(blast_results)
}

empty_blast_tsv_tab <- function() {
  tab <- dplyr::data_frame(path = character(0), query = character(0), gi = character(0))
  class(tab) <- append(class(tab), "blast_tabular")
  return(tab)
}

parse_blast_tsv_safe <- purrr::safely(parse_blast_tsv, otherwise = empty_blast_tsv_tab())

gi2taxid <- function(tab, taxdb) {
  UseMethod("gi2taxid", tab)
}

#' @param tab data_frame with parsed blast results
#' @param taxdb sqlite database with gi_taxid_nucl and gi_taxid_prot tables
gi2taxid.blast_tabular <- function(tab, taxdb) {

  mapped_gis <- tab$gi

  message("Connect to database")
  db <- dplyr::src_sqlite(taxdb, create = TRUE)

  message("Collect tax_ids from tables")
  gi_nuc <- dplyr::tbl(db, "gi_taxid_nucl") %>%
    dplyr::filter(gi %in% mapped_gis) %>%
    dplyr::collect()
  gi_prot <- dplyr::tbl(db, "gi_taxid_prot") %>%
    dplyr::filter(gi %in% mapped_gis) %>%
    dplyr::collect()

  message("Bind rows")
  gi_tab <- dplyr::bind_rows(gi_nuc, gi_prot) %>%
    dplyr::mutate_at(dplyr::vars(gi), as.character)

  message("Join taxonomy to blast results by gi")
  known <- dplyr::left_join(tab, gi_tab)

  message("Fill in few missing tax_ids by quering remote ncbi database")
  with_taxid <- dplyr::filter(known, !is.na(tax_id))
  no_taxid <- dplyr::filter(known, is.na(tax_id))

  query <- dplyr::mutate(no_taxid, tax_id = map_chr(gi, query_taxid)) %>%
    dplyr::select(gi, tax_id) %>%
    dplyr::mutate_at("tax_id", as.integer)

  fixed_taxid <- dplyr::full_join(dplyr::select(no_taxid, -tax_id), query)
  blast_results_taxids <- dplyr::bind_rows(with_taxid, fixed_taxid)
  class(blast_results_taxids) <- append(class(blast_results_taxids), "blast_results_taxids")
  return(blast_results_taxids)
}

filter_division <- function(tab, nodes, div_id, div, not_div) {
  UseMethod("filter_division", tab)
}

#' @param tab blast results tab with tax_ids, a data_frame.
#' @param nodes path to nodes.csv file, a character string.

#' @param div path ot output.csv file for records belonging to div_id, a character string.
#' @param not_div path ot output.csv file for records NOT belongigng to div_id, a character string.
filter_division.blast_results_taxids <- function(tab, nodes, div_id, div, not_div) {

  message("Import gi tax_id table and taxonomy nodes table")
  nodes <- readr::read_csv(nodes)

  message("Merge names by tax_id to known sequences")
  mapped_tab <- dplyr::left_join(tab, nodes)

  division <- dplyr::filter(mapped_tab, division_id == div_id)
  not_division <- dplyr::filter(mapped_tab, division_id != div_id)

  readr::write_csv(division, div)
  readr::write_csv(not_division, not_div)
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

  blast_results <- grep("tsv$", dots, value = TRUE)
  tab <- parse_blast_tsv_safe(blast_results)

  # if it's not empty result then stop error
  if (!is.null(tab$error) && !stringr::str_detect(tab$error$message, "No BLAST hits to parse")) stop(tab$error)

  message("Map tax_ids to gis\n")
  tab <- gi2taxid(tab$result, taxdb)

  message("dplyr::filter phages (division_id == 3) and save results to csv files\n")
  filter_division(tab = tab, nodes = nodes, div = phages, not_div = viruses, div_id = div_id)
}

do.call(blast_taxonomy, c(snakemake@input, snakemake@output))
