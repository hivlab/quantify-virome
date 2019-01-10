
library(magrittr)

#' Query taxonomy id number for nucleotide GenInfo Identifier (GI) from NCBI nucleotide database
#' @param gi GenInfo Identifier (GI), a character string.
query_taxid <- function(gi) {
  res <- httr::GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
             query = list(db = "nucleotide", id = gi, rettype = "fasta", retmode = "xml"))
  cont <- httr::content(res, as = "parsed", encoding = "UTF-8")
  rvest::xml_node(cont, xpath = "//TSeq_taxid") %>% xml2::xml_text()
}

#' @param gids GenInfo Identifier (GI), a character (or integer) vector.
get_taxid <- function(gids) {
  if (Sys.getenv("NCBI_API_KEY") == "") {
    warning("NCBI_API_KEY environment variable not set.\nPosting more than 3 requests per second to the NCBI E-utilities without an API key will receive an error message. API key can be obtained from the Settings page of their NCBI account (to create an account, visit http://www.ncbi.nlm.nih.gov/account/).\n")
    api_key <- NULL
  } else {
    api_key <- Sys.getenv("NCBI_API_KEY")
  }
  res <- httr::GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", query = list(db = "nucleotide", id = paste0(gids, collapse = ";"), rettype = "fasta", retmode = "xml", api_key = api_key))
  cont <- httr::content(res, as = "parsed", encoding = "UTF-8")
  xml2::xml_children(cont) %>%
    purrr::map(xml2::xml_find_first, ".//TSeq_taxid") %>%
    purrr::map(xml2::xml_text) %>%
    unlist()
}

#' Import and parse tabular BLAST+ output (outfmt 6) to tibble
#' @param ... path(s) to tabular BLAST+ output(s), a character vector.
parse_blast_tsv <- function(...) {

  message("Importing BLAST+ tabular output")
  blast_results <- dplyr::tibble(path = c(...)) %>%
    dplyr::mutate(tsv = purrr::map(path, readr::read_tsv)) %>%
    tidyr::unnest()

  message("Munging metadata")
  blast_results <- dplyr::rename(blast_results, query = qseqid, gi = sgi) %>%
    dplyr::mutate_at(dplyr::vars(gi), as.character) %>%
    dplyr::mutate(path = basename(path)) %>%
    dplyr::select(path, query, gi, dplyr::everything())
  class(blast_results) <- append(class(blast_results), "blast_tabular")
  return(blast_results)
}

#' Create empty data_frame of class blast_tabular with path and all std outfmt 6 columns
empty_blast_tsv_tab <- function() {
  tab <- dplyr::tibble(path = character(0), query = character(0), gi = character(0), pident = numeric(0),
                           length = integer(0), mismatch = integer(0), gapopen = integer(0), qstart = integer(0),
                           qend = integer(0), sstart = integer(0), send = integer(0),
                           evalue = numeric(0), bitscore = numeric(0))
  class(tab) <- append(class(tab), "blast_tabular")
  return(tab)
}

#' Parse blast output safely, in case of error output empty blast_tabular dataframe
parse_blast_tsv_safe <- purrr::safely(parse_blast_tsv, otherwise = empty_blast_tsv_tab())

gi2taxid <- function(tab, taxdb, api_key) {
  UseMethod("gi2taxid", tab)
}

#' @param tab data_frame with parsed blast results of class blast_tabular.
#' @param taxdb sqlite database with gi_taxid_nucl and gi_taxid_prot tables.
#' @param api_key ncbi API key, a character string, defaults to NULL.
gi2taxid.blast_tabular <- function(tab, taxdb, api_key) {

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

  query <- dplyr::mutate(no_taxid, tax_id = get_taxid(gi)) %>%
    dplyr::select(gi, tax_id) %>%
    dplyr::mutate_at("tax_id", as.integer)

  fixed_taxid <- dplyr::full_join(dplyr::select(no_taxid, -tax_id), query)
  blast_results_taxids <- dplyr::bind_rows(with_taxid, fixed_taxid)
  class(blast_results_taxids) <- append(class(blast_results_taxids), "blast_results_taxids")
  return(blast_results_taxids)
}

filter_division <- function(tab, nodes, division_id, div, not_div) {
  UseMethod("filter_division", tab)
}

#' Filter BLAST hits belonging to div_id
#' @param tab blast results tab with tax_ids, a data_frame of class blast_results_taxids.
#' @param nodes path to nodes.csv file, a character string.
#' @param division_id taxonomic division id, an integer vector.
#' @param div path ot output.csv file for records belonging to div_id, a character string.
#' @param not_div path ot output.csv file for records NOT belongigng to div_id, a character string.
filter_division.blast_results_taxids <- function(tab, nodes, division_id, div, not_div) {

  message("Import gi tax_id table and taxonomy nodes table")
  nodes <- readr::read_csv(nodes)

  message("Merge names by tax_id to known sequences")
  mapped_tab <- dplyr::left_join(tab, nodes)

  # Just in case, ensure that we have division_ids as integers
  div_id <- as.integer(unlist(division_id))

  # Filter hits by division_id
  is_division <- mapped_tab %>%
    dplyr::group_by(query) %>%
    dplyr::filter(all(division_id %in% div_id))
  not_division <- dplyr::anti_join(mapped_tab, is_division)

  # Write results to files
  readr::write_csv(is_division, div)
  readr::write_csv(not_division, not_div)
}

#' Parse and add taxonomy to blast+ results
#' @param ... paths(s) to tabular blast results file (outfmt 6), a chracter string.
#' @param taxdb path to taxonomy database, sqlite database, a character string.
#' @param nodes path to taxonomy nodes.csv file, a chracter string.
#' @param division path to output file for hits assigned to division(s) of interest, a character string.
#' @param other path to output file for other hits not belonging to division(s) of interest, character string
#' @param division_id division id of interest, integer vector. Defaults to 3, phages.
#' @param api_key ncbi API key, a character string, defaults to NULL.
#'
blast_taxonomy <- function(..., taxdb, nodes, division, other, division_id = 3, api_key = NULL) {

  message("Parse blast results\n")
  dots <- c(...)

  blast_results <- grep("tsv$", dots, value = TRUE)
  tab <- parse_blast_tsv_safe(blast_results)

  # If it's not empty result then stop with error
  if (!is.null(tab$error) && !stringr::str_detect(tab$error$message, "object 'qseqid' not found")) stop(tab$error)

  message("Map tax_ids to gis\n")
  tab <- gi2taxid(tab$result, taxdb, api_key)

  message("Filter phages (division_id == 3) and save results to csv files\n")
  filter_division(tab = tab, nodes = nodes, division_id = division_id, div = division, not_div = other)
}

do.call(blast_taxonomy, c(snakemake@input, snakemake@output, snakemake@params))
