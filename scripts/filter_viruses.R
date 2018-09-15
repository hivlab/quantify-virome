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

  data_frame(query, gi, Hit_accession, Hit_def) %>% bind_cols(hsps)
}

query_taxid <- function(gi) {
  res <- GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
             query = list(db = "nucleotide", id = gi, rettype = "fasta", retmode = "xml"))
  cont <- content(res, as = "parsed", encoding = "UTF-8")
  xml_node(cont, xpath = "//TSeq_taxid") %>% xml_text()
}

blast_taxonomy <- function(nt_virus, nr_virus, taxdb, nodes, phages, viruses) {

  message("Open and fix BLAST+ xml strings")
  virus <- data_frame(path = c(nt_virus, nr_virus)) %>%
    mutate(xml = map(path, read_lines)) %>%
    mutate(xml = map(xml, split_hits)) %>%
    unnest() %>%
    mutate(xml = map(xml, str_c, collapse = ""))


  message("Import fixed xml")
  virus <- mutate(virus, xml = map(xml, read_xml))

  message("Extract hits from xml")
  virus <- mutate(virus, hits = map(xml, parse_hits))

  message("Extract blast results from list")
  virus <- virus %>%
    dplyr::select(path, hits) %>%
    unnest() %>%
    mutate(blast = case_when(
      str_detect(path, "blastn") ~ "blastn",
      str_detect(path, "blastx") ~ "blastx"
    )) %>%
    dplyr::select(-path)

  message("Query local vhunter database")
  db <- src_sqlite(taxdb, create = TRUE)
  message("tbl")
  nucl_db <- tbl(db, "gi_taxid_nucl")
  prot_db <- tbl(db, "gi_taxid_prot")
  message("collect")
  gi_nuc <- nucl_db %>%
    filter(gi %in% virus$gi) %>%
    collect()
  gi_prot <- prot_db %>%
    filter(gi %in% virus$gi) %>%
    collect()
  message("bind rows")
  gi_tab <- bind_rows(gi_nuc, gi_prot) %>%
    mutate_at(vars(gi), as.character)

  message("Join taxonomy by gi")
  known <- left_join(virus, gi_tab)

  message("Fill in few missing tax_ids by quering remote ncbi database")
  with_taxid <- filter(known, !is.na(tax_id))
  no_taxid <- filter(known, is.na(tax_id))

  query <- mutate(no_taxid, tax_id = map_chr(gi, query_taxid)) %>%
    dplyr::select(gi, tax_id) %>%
    mutate_at("tax_id", as.integer)

  fixed_taxid <- full_join(dplyr::select(no_taxid, -tax_id), query)
  known_fix <- bind_rows(with_taxid, fixed_taxid)

  message("Import gi tax_id table and taxonomy nodes table")
  nodes <- read_csv(nodes)

  message("Merge names by tax_id to known sequences")
  known_fix <- left_join(known_fix, nodes)

  phage <- filter(known_fix, near(division_id, 3))
  candidate_viruses <- filter(known_fix, !near(division_id, 3))

  write_csv(phage, phages)
  write_csv(candidate_viruses, viruses)
}

blast_taxonomy(nt_virus = snakemake@input[[1]],
               nr_virus = snakemake@input[[2]],
               taxdb = snakemake@input[[3]],
               nodes = snakemake@input[[4]],
               phages = snakemake@output[[1]],
               viruses = snakemake@output[[2]])
