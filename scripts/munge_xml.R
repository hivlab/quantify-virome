library(tidyverse)
library(xml2)

# Blast output against virus database
nt_virus <- "output/I1164-mini_blastn_virus_1_known-viral.xml"
nr_virus <- "output/I1164-mini_blastx_virus_1_known-viral.xml"

# Fix BLAST+ xml
split_hits <- function(xml) {
  n <- diff(c(which(str_detect(xml, "<\\?xml")), length(xml) + 1))
  split(xml, rep(1:length(n), n))
}

virus <- data_frame(path = c(nt_virus, nr_virus)) %>%
  mutate(xml = map(path, read_lines))
virus <- mutate(virus, xml = map(xml, split_hits))
virus <- mutate(virus, xml = map(xml, ~map(.x, str_c, collapse = "")))

# Import fixed xml
virus <- mutate(virus, xml = map(xml, map, read_xml))

# Convert xml to list
virus <- mutate(virus, xml_list = map(xml, map, as_list))

# Extract blast results from list
xml_hits <- virus %>%
  unnest() %>%
  mutate(BlastOutput = map(xml_list, "BlastOutput"),
         iteration = map(BlastOutput, "BlastOutput_iterations"),
         iteration = map(iteration, "Iteration"),
         query = map_chr(iteration, ~.x$"Iteration_query-def"[[1]]),
         hits = map(iteration, "Iteration_hits")) %>%
  select(-xml, -xml_list)

parse_hits <- function(hits) {
  ids <-  map_dfr(hits, `[`, c("Hit_id", "Hit_accession", "Hit_def")) %>%
    unnest() %>%
    mutate(gi = str_replace(Hit_id, "gi\\|(\\d+)\\|.*", "\\1")) %>%
    select(gi, Hit_accession, Hit_def)
  hsps <- map(hits, ~.x$Hit_hsps$Hsp) %>%
    map_dfr(`[`, c("Hsp_evalue", "Hsp_bit-score", "Hsp_score", "Hsp_query-from",
                       "Hsp_hit-from", "Hsp_hit-to", "Hsp_query-frame", "Hsp_hit-frame",
                       "Hsp_query-to","Hsp_identity", "Hsp_positive", "Hsp_align-len",
                       "Hsp_gaps", "Hsp_qseq", "Hsp_hseq", "Hsp_midline")) %>%
    unnest()
  bind_cols(ids, hsps)
}

xml_hits <- mutate(xml_hits, hits = map(hits, parse_hits))

xml_hits_unnested <- xml_hits %>%
  select(path, query, hits) %>%
  unnest() %>%
  mutate(blast = case_when(
    str_detect(path, "blastn") ~ "blastn",
    str_detect(path, "blastx") ~ "blastx"
  ))

# Use vhunter database
db_path <- "/gpfs/software/VirusSeeker/databases/taxdump_300817/vhunter.db"
db <- src_sqlite(db_path, create = TRUE)
nucl_db <- tbl(db, "gi_taxid_nucl")
prot_db <- tbl(db, "gi_taxid_prot")
gi_nuc <- nucl_db %>%
  filter(gi %in% xml_hits_unnested$gi) %>%
  collect()
gi_prot <- prot_db %>%
  filter(gi %in% xml_hits_unnested$gi) %>%
  collect()
gi_tab <- bind_rows(gi_nuc, gi_prot)

# Join taxonomy by gi
known_unnested <- left_join(mutate_at(xml_hits_unnested, "gi", as.integer), gi_tab)
no_taxid <- filter(known_unnested, is.na(tax_id))

library(httr)
library(xml2)
library(rvest)
query_taxid <- function(gi) {
  res <- GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
             query = list(db = "nucleotide", id = gi, rettype = "fasta", retmode = "xml"))
  cont <- content(res, as = "parsed")
  xml_node(cont, "TSeq_taxid") %>% xml_text()
}

query <- mutate(no_taxid, tax_id = map_chr(gi, query_taxid)) %>%
  select(gi, tax_id) %>%
  mutate_at("tax_id", as.integer)

 <- left_join(no_taxid, query)

filter(no_taxid, is.na(tax_id))

